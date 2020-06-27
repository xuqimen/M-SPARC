function [upper_bound_guess_vecs,psi,EigVal,a0,bup,lambda_cutoff,ChebComp] = eigSolver_sq(S,count,upper_bound_guess_vecs,psi,EigVal,a0,bup,lambda_cutoff)
% @brief    Solve the linearized Hamiltonian eigenproblem.
%
% @param count                  : SCF iteration count.
% @param upper_bound_guess_vecs : guess vector for maximum eigenvector of
%                                 the linearized Hamiltonian.
% @param psi                    : eigenvectors of the previous linearized
%                                 Hamiltonian.
% @param EigVal                 : previous eigenvalues
% @param a0                     : lower bound for Chebyshev filtering (1st SCF).
% @param bup                    : upper bound for Chebyshev filtering (1st SCF).
% @param lambda_cutoff          : cut-off for chebyshev filtering (1st SCF).
%
% @authors	Qimen Xu <qimenxu@gatech.edu>
%           Abhiraj Sharma <asharma424@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech
%

if S.ForceCount > 1
	S.rhoTrigger = 1;
end

if S.parallel ~= 1
    for ks = 1:S.tnkpt*S.nspin
		if ks <= S.tnkpt
			kpt = ks;
			spin = 1;
		else
			kpt = ks - S.tnkpt;
			spin = 2;
		end
		% Heff = spdiags(S.Veff(:,spin),0,S.N,S.N);
		Heff = S.Veff(:,spin);
        rng('default'); % Initialize random number generator
		rng(ks+1);
		%opts = struct('maxit', 10000, 'tol', 1e-6, 'p', S.Nev+10, 'v0', rand(S.N,1), 'isreal', true);
		opts = struct('maxit', 100, 'tol', S.TOL_LANCZOS, 'v0', rand(S.N,1));
		kpt_vec = S.kptgrid(kpt,:);
		[DL11,DL22,DL33,DG1,DG2,DG3] = blochLaplacian_1d(S,kpt_vec);
		Hfun = @(x) h_nonlocal_vector_mult(DL11,DL22,DL33,DG1,DG2,DG3,Heff,x,S,kpt_vec);
		if ~(isreal(DL11) && isreal(DL22) && isreal(DL33))
			opts.isreal = false;
		end

		if(count == 1)
			% Upper bound estimator
			opts.maxit = 300; % WARNING: might need more accuracy
			if(S.ForceCount == 1)
				% For first relaxation step
				[upper_bound_guess_vecs(:,ks), bup(ks)] = (eigs(Hfun,S.N,1,'lr',opts));
				bup(ks) = real(bup(ks)) * 1.01;
			else
				% For subsequent relaxation steps
				opts.v0 = S.upper_bound_guess_vecs(:,ks);
				[upper_bound_guess_vecs(:,ks), bup(ks)] = eigs(Hfun,S.N,1,'lr',opts) ;
				bup(ks) = real(bup(ks)) * 1.01;
			end
			% Lower bound estimator
			if(S.ForceCount == 1)
				% For first relaxation step
				a0(ks) = real(eigs(Hfun,S.N,1,'sr',opts)) - 0.1;
			else
				% For subsequent relaxation steps use lowest eigenvalue
				a0(ks) = S.EigVal(1,ks);
			end
			% Always use this on every first SCF of relaxation
			%lambda_cutoff(ks) = 0.5 * (bup(ks) + a0(ks));
		else
			% For subsequent steps
			if (count > S.rhoTrigger)
				if S.chefsibound_flag == 1
					% Upper bound
					opts.tol = S.TOL_LANCZOS; % WARNING: might need more accuracy than the default
					opts.v0 = S.upper_bound_guess_vecs(:,ks);
					[upper_bound_guess_vecs(:,ks), bup(ks)] = eigs(Hfun,S.N,1,'lr',opts);
					bup(ks) = real(bup(ks)) * 1.01;
				end
				% Lower bound
				a0(ks) = S.EigVal(1,ks);
			end
			% Set the filter cut off
			%S.lambda_f + log10(1e6-1)/S.bet + 0.5; 
			%lambda_cutoff(ks) = S.EigVal(S.Nev,ks) + 0.10; 
		end

		if count == 1 && S.ForceCount == 1
			lambda_cutoff(ks) = 0.5 * (bup(ks) + a0(ks));
		else
			% Set the filter cut off
			%S.lambda_f + log10(1e6-1)/S.bet + 0.5; 
			lambda_cutoff(ks) = S.EigVal(S.Nev,ks) + 0.10; 
		end
		
		fprintf('filter cutoff = %f, lower bound = %f, upper bound = %f\n',lambda_cutoff(ks),a0(ks),bup(ks));
		% Chebyshev filtering
		psi(:,:,ks) = chebyshev_filter(psi(:,:,ks),S.npl,lambda_cutoff(ks),bup(ks),a0(ks),DL11,DL22,DL33,DG1,DG2,DG3,Heff,S,kpt_vec);
		%psi(:,:,ks) = orth(psi(:,:,ks));
		psi(:,:,ks) = OrthChol(psi(:,:,ks));
		Nev1 = size(psi(:,:,ks),2);    % WARNING: ORTH(psi) might change the size of psi, that's why we update Nev
		assert(Nev1 == S.Nev,'Number of states have changed within SCF');

		% Subspace Hamiltonian
		Hs = psi(:,:,ks)' * h_nonlocal_vector_mult(DL11,DL22,DL33,DG1,DG2,DG3,Heff,psi(:,:,ks),S,kpt_vec);

		% Solve subspace eigenproblem,
		if S.cell_typ < 3
			Hs = 0.5 * (Hs + Hs');
		end
		
		if S.densMatFlag == 1
			% instead of solving the eigenproblem, we find the subspace density
			% matrix from Hs: Ds = f(Hs,Ef,sigma), since Ef is unknown, we
			% first find the Tj((Hs-c)/e) matrices and save them (just need
			% to save the trace for the calculation of Ef later, and once
			% Ef is known, recalculate Tj's again, but may be slower)
			opts_Hs.maxit = 300;
			opts_Hs.tol = 1e-2;
			opts_Hs.disp = 1;
			tic
% 			EigVal(1    ,ks) = real(eigs(Hs,1,'sr',opts_Hs));
% 			EigVal(S.Nev,ks) = real(eigs(Hs,1,'lr',opts_Hs));
			[Lmin,Lmax] = Lanczos_Algo(Hs,1e-10);
			EigVal(1,ks) = Lmin;
			EigVal(S.Nev,ks) = Lmax;
			fprintf('eigmin(Hs) = %f, eigmax(Hs) = %f\n',EigVal(1,ks),EigVal(S.Nev,ks));
			toc
			
			tic
			eigmin = EigVal(1,ks)-0.2-0.1*abs( EigVal(1,ks));
			eigmax = EigVal(S.Nev,ks)+0.2+0.1*abs( EigVal(S.Nev,ks));
			ChebComp = Chebyshev_matvec_comp(S.sq_npl,Hs,eye(Nev1),eigmin,eigmax);
			ChebComp(1).eigmin(ks) = eigmin; % needs to be exactly the same as the one used in Chebyshev_matvec_comp
			ChebComp(1).eigmax(ks) = eigmax; % needs to be exactly the same as the one used in Chebyshev_matvec_comp
			ChebComp(1).sq_npl = S.sq_npl;
			toc
		else
			[Q, Q1] = eig(Hs);
			EigVal(:,ks) = real(diag(Q1)); % WARNING: Taking real part only!
			% subspace rotation
			psi(:,:,ks) = psi(:,:,ks) * Q;
			% Normalize psi, s.t. integral(psi_new' * psi_new) = 1
			scfac = 1 ./ sqrt(sum(repmat(S.W,1,S.Nev) .* (psi(:,:,ks) .* conj(psi(:,:,ks))),1));
			% psi(:,:,ks) = psi(:,:,ks) * diag(scfac);
			psi(:,:,ks) = bsxfun(@times, psi(:,:,ks), scfac);
		end
    end
else 
	% Before getting into parfor, set to use only one thread
	LASTN = maxNumCompThreads(1);
	parfor (ks = 1:S.tnkpt*S.nspin, S.num_worker_heuristic)
		if ks <= S.tnkpt
			kpt = ks;
			spin = 1;
		else
			kpt = ks - S.tnkpt;
			spin = 2;
		end
		% Heff = spdiags(S.Veff(:,spin),0,S.N,S.N);
		Heff = S.Veff(:,spin);
        rng('default'); % Initialize random number generator
		rng(ks+1);
		%opts = struct('maxit', 10000, 'tol', 1e-6, 'p', S.Nev+10, 'v0', rand(S.N,1), 'isreal', true);
		opts = struct('maxit', 100, 'tol', S.TOL_LANCZOS, 'v0', rand(S.N,1));
		kpt_vec = S.kptgrid(kpt,:);
		[DL11,DL22,DL33,DG1,DG2,DG3] = blochLaplacian_1d(S,kpt_vec);
		Hfun = @(x) h_nonlocal_vector_mult(DL11,DL22,DL33,DG1,DG2,DG3,Heff,x,S,kpt_vec);
		if ~(isreal(DL11) && isreal(DL22) && isreal(DL33))
			opts.isreal = false;
		end

		if (count == 1)
			% Upper bound estimator
			opts.maxit = 300; % WARNING: might need more accuracy
			if(S.ForceCount == 1)
				% For first relaxation step
				[upper_bound_guess_vecs(:,ks), bup(ks)] = (eigs(Hfun,S.N,1,'lr',opts));
				bup(ks) = real(bup(ks)) * 1.01;
			else
				% For subsequent relaxation steps
				opts.v0 = S.upper_bound_guess_vecs(:,ks);
				[upper_bound_guess_vecs(:,ks), bup(ks)] = eigs(Hfun,S.N,1,'lr',opts) ;
				bup(ks) = real(bup(ks)) * 1.01;
			end
			% Lower bound estimator
			if(S.ForceCount == 1)
				% For first relaxation step
				a0(ks) = real(eigs(Hfun,S.N,1,'sr',opts)) - 0.1;
			else
				% For subsequent relaxation steps use lowest eigenvalue
				a0(ks) = S.EigVal(1,ks);
			end
			% Always use this on every first SCF of relaxation
			%lambda_cutoff(ks) = 0.5 * (bup(ks) + a0(ks));
		else
			% For subsequent steps
			if (count > S.rhoTrigger)
				if S.chefsibound_flag == 1
					% Upper bound
					opts.tol = S.TOL_LANCZOS; % WARNING: might need more accuracy than the default
					opts.v0 = S.upper_bound_guess_vecs(:,ks);
					[upper_bound_guess_vecs(:,ks), bup(ks)] = eigs(Hfun,S.N,1,'lr',opts);
					bup(ks) = real(bup(ks)) * 1.01;
				end
				% Lower bound
				a0(ks) = S.EigVal(1,ks);
			end
			% Set the filter cut off
			%S.lambda_f + log10(1e6-1)/S.bet + 0.5; 
			%lambda_cutoff(ks) = S.EigVal(S.Nev,ks) + 0.10; 
		end
		
		if count == 1 && S.ForceCount == 1
			lambda_cutoff(ks) = 0.5 * (bup(ks) + a0(ks));
		else
			% Set the filter cut off
			%S.lambda_f + log10(1e6-1)/S.bet + 0.5; 
			lambda_cutoff(ks) = S.EigVal(S.Nev,ks) + 0.10; 
		end
		
		% fprintf('filter cutoff = %f, lower bound = %f, upper bound = %f\n',lambda_cutoff(ks),a0(ks),bup(ks));
		% Chebyshev filtering
		psi(:,:,ks) = chebyshev_filter(psi(:,:,ks),S.npl,lambda_cutoff(ks),bup(ks),a0(ks),DL11,DL22,DL33,DG1,DG2,DG3,Heff,S,kpt_vec);
		psi(:,:,ks) = orth(psi(:,:,ks));
		Nev1 = size(psi(:,:,ks),2);    % WARNING: ORTH(psi) might change the size of psi, that's why we update Nev
		assert(Nev1 == S.Nev,'Number of states have changed within SCF');

		% Subspace Hamiltonian
		Hs = psi(:,:,ks)' * h_nonlocal_vector_mult(DL11,DL22,DL33,DG1,DG2,DG3,Heff,psi(:,:,ks),S,kpt_vec);

		% Solve subspace eigenproblem,
		if S.cell_typ < 3
			Hs = 0.5 * (Hs + Hs');
		end
		[Q, Q1] = eig(Hs);
		EigVal(:,ks) = real(diag(Q1)); % WARNING: Taking real part only!

		% subspace rotation
		psi(:,:,ks) = psi(:,:,ks) * Q;

		% Normalize psi, s.t. integral(psi_new' * psi_new) = 1
		scfac = 1 ./ sqrt(sum(repmat(S.W,1,S.Nev) .* (psi(:,:,ks) .* conj(psi(:,:,ks))),1));
		% psi(:,:,ks) = psi(:,:,ks) * diag(scfac);
		psi(:,:,ks) = bsxfun(@times, psi(:,:,ks), scfac);
    end

	% After parfor, reset the number of threads as before
	maxNumCompThreads(LASTN);
end

end


function ChebComp = Chebyshev_matvec_comp(npl,A,X,a,b)
% Chebyshev_matvec_comp evaluates Ti(A,a,b) for i = 0,...,npl, where Ti(x,a,b)
% is the standard Chebyshev polynomial Ti(x) with x in [a,b] mapped to [-1,1].
% The results Ti(A,a,b) are saved in a structure ChebComp(i).Ti
a0 = a;
e = (b-a)/2;
c = (b+a)/2;
%sigma = e/(a0 - c); 
sigma = e/(c - a0); 
sigma1 = sigma;
gamma = 2/sigma1;

ChebComp(1).Hs = A; % TODO: REMOVE AFTER CHECK

%Ysum = C(1) * X;
ChebComp(1).Ti = X;
ChebComp(1).tr_Ti = trace(ChebComp(1).Ti);
if (npl <= 0), return; end
%A(1:size(A,1)+1:end) = A(1:size(A,1)+1:end) - c;
%Y = (sigma1/e)*(A*X); % T1(A,a,b) * X
Y = (sigma1/e)*(A*X-c*X); % T1(A,a,b) * X
ChebComp(2).Ti = Y;
ChebComp(2).tr_Ti = trace(ChebComp(2).Ti);
ii = 2;
while(ii <= npl)
	sigma2 = 1/(gamma - sigma);
	% AX = A * Y;
	Ynew = (2*sigma2/e)*(A*Y - c*Y) - ((sigma*sigma2)*X);
	%Ynew = (2*sigma2/e)*(A*Y) - ((sigma*sigma2)*X);
	%Ysum = Ysum + C(ii+1) * Ynew;
	ChebComp(ii+1).Ti = Ynew;
	ChebComp(ii+1).tr_Ti = trace(ChebComp(ii+1).Ti);
	X = Y;
	Y = Ynew;
	sigma = sigma2;
	ii = ii + 1;
end

end