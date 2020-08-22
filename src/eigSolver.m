function [upper_bound_guess_vecs,psi,EigVal,a0,bup,lambda_cutoff,CS] = eigSolver(S,count,upper_bound_guess_vecs,psi,EigVal,a0,bup,lambda_cutoff)
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

% save Complemenary subspace related variables in CS struct
if S.CSFlag == 1
	CS.psi_t = zeros(S.N,S.Ns_top,S.tnkpt*S.nspin);
else 
	CS = struct([]);
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
				a0(ks) = min(S.EigVal(:,ks));
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
				a0(ks) = min(S.EigVal(:,ks));
			end
			% Set the filter cut off
			%S.lambda_f + log10(1e6-1)/S.bet + 0.5; 
			%lambda_cutoff(ks) = max(S.EigVal(:,ks)) + 0.10; 
		end

		if count == 1 && S.ForceCount == 1
			lambda_cutoff(ks) = 0.5 * (bup(ks) + a0(ks));
		else
			% Set the filter cut off
			%S.lambda_f + log10(1e6-1)/S.bet + 0.5; 
			lambda_cutoff(ks) = max(S.EigVal(:,ks)) + 0.10; 
		end
		
		% fprintf('filter cutoff = %f, lower bound = %f, upper bound = %f\n',lambda_cutoff(ks),a0(ks),bup(ks));
		% Chebyshev filtering
		psi(:,:,ks) = chebyshev_filter(psi(:,:,ks),S.npl,lambda_cutoff(ks),bup(ks),a0(ks),DL11,DL22,DL33,DG1,DG2,DG3,Heff,S,kpt_vec);
		%psi(:,:,ks) = orth(psi(:,:,ks)); % uses svd
		psi(:,:,ks) = orthChol(psi(:,:,ks)); % uses Cholesky
		Nev1 = size(psi(:,:,ks),2);    % WARNING: ORTH(psi) might change the size of psi, that's why we update Nev
		assert(Nev1 == S.Nev,'Number of states have changed within SCF');

		% Subspace Hamiltonian
		Hs = psi(:,:,ks)' * h_nonlocal_vector_mult(DL11,DL22,DL33,DG1,DG2,DG3,Heff,psi(:,:,ks),S,kpt_vec);

		% Solve subspace eigenproblem,
		if S.cell_typ < 3
			Hs = 0.5 * (Hs + Hs');
		end
		
		% S.CSFlag = 1;
		% S.Ns_top = round(S.Nev * 0.2);
		% S.Ns_top = S.Nev * 2 - S.Nelectron;
		
		if S.CSFlag == 1 % use Complementary Subspace method to solve the subspace eigenproblem
			CS.tr_Hs = trace(Hs);
% 			CS_index_reverse = S.CS_index(end:-1:1); 
			t1 = tic;
			
			
			%a_CS = a0(ks) - 0.2;
			%a0_CS = lambda_cutoff(ks);
			
			%------------ Test worked -------------------
% 			a_CS = lambda_cutoff(ks);
% 			a0_CS = a0(ks) - 0.2;
% 			if count == 1 && S.ForceCount == 1
% 				b_CS = (a_CS + a0_CS)*0.5;
% 				S.Qt = rand(S.Nev,S.Ns_top) - 0.5;
% 			else
% 				b_CS = EigVal(S.Nev-S.Ns_top+1,ks) - 0.20;
% 			end
			
% 			%------------ Test 2 -------------------
% 			a0_CS = -lambda_cutoff(ks);
% 			b_CS = -(a0(ks) - 0.2);
% 			if count == 1 && S.ForceCount == 1
% 				a_CS = (b_CS + a0_CS)*0.5;
% 				S.Qt = rand(S.Nev,S.Ns_top) - 0.5;
% 			else
% 				a_CS = -(EigVal(S.Nev-S.Ns_top+1,ks) - 0.20);
% 			end

			%------------ Test 3 -------------------
% 			a0_CS = lambda_cutoff(ks);
% 			b_CS = (a0(ks) - 0.2);
% 			if count == 1 && S.ForceCount == 1
% 				a_CS = (b_CS + a0_CS)*0.5;
% 				S.Qt = rand(S.Nev,S.Ns_top) - 0.5;
% 			else
% 				a_CS = (EigVal(S.Nev-S.Ns_top+1,ks) - 0.20);
% 			end
			
			%------------ Test 4 -------------------
% 			a0_CS = -(a0(ks) - 0.2);
% 			b_CS = -lambda_cutoff(ks);
% 			if count == 1 && S.ForceCount == 1
% 				a_CS = -(b_CS + a0_CS)*0.5;
% 				S.Qt = rand(S.Nev,S.Ns_top) - 0.5;
% 			else
% 				a_CS = -(EigVal(S.Nev-S.Ns_top+1,ks) - 0.20);
% 			end
			
			%------------ Test 5 -------------------
			a0_CS = lambda_cutoff(ks)-0.1;
			b_CS = (a0(ks) - 0.05);
			if count == 1 && S.ForceCount == 1
				a_CS = (b_CS + a0_CS)*0.5;
				S.Qt = rand(S.Nev,S.Ns_top) - 0.5;
			else
				a_CS = (EigVal(S.Nev-S.Ns_top+1,ks) - 0.05);
			end
			m = 16; n = 1;
			[Qt_temp, Dt_temp] = eigsolver_ChebFilt(Hs, S.Qt, m, n, a_CS, b_CS, a0_CS);
			Qt = Qt_temp;
			Dt = Dt_temp;
			disp('Using 2nd level of CheFSI')
			diag(Dt_temp)
			
			% Method 2: directly use eigs to find top Nt states
			% [Qt, Dt] = eigs(Hs,S.Ns_top,'lr');
			% [eigval_sort,ind_sort] = sort(diag(Dt)); % sort eigvals to be increasingly ordered
			% Qt = Qt(:,ind_sort);
			% Dt = diag(eigval_sort);

			CS.Qt = Qt_temp; % save top eigenvectors of Hs for next iteration
			CS.Hs = Hs; % TODO: remove after check
			fprintf('Time for solving top %d eigenvalues of Hs: %f s\n',S.Ns_top,toc(t1));

			
			%EigVal(CS_index_reverse,ks) = real(diag(Dt)); % since here eigenvalues are descendingly ordered
			EigVal(S.CS_index,ks) = sort(real(diag(Dt))); % since here eigenvalues are descendingly ordered
			
			% find minimum eigenvalue as a CheFSI bound
			[Lmin,~] = lanczos_algo(Hs,1e-2);
			EigVal(1,ks) = Lmin;
			fprintf('eigmin(Hs) = %f, eigmax(Hs) = %f\n',EigVal(1,ks),EigVal(S.Nev,ks));
			
			EigVal(2:S.Nev-S.Ns_top) = nan; % TODO: remove after check
			EigVal(:,ks)
			
			Q1_temp = eig(Hs);
			disp('Actual eigvals:')
			disp(Q1_temp)
			fprintf('Sum of eigvals = %.15f, Trace of Hs = %.15f, diff = %.3e\n',...
				sum(Q1_temp), CS.tr_Hs, abs(sum(Q1_temp)-CS.tr_Hs));
			
			size(Q1_temp)
			size(diag(Dt_temp))
			fprintf('Accuracy of 2nd CheFSI eigensolver: %.3e\n',abs(Q1_temp(S.CS_index) - diag(Dt_temp)));
			
			
			% rotate the top states and save result in psi_t (NOTE: may need to update psi for next SCF)
			%psi(:,CS_index_reverse,ks) = psi(:,:,ks) * Qt;
			%CS.psi_t(:,:,ks) = psi(:,:,ks) * Qt(:,end:-1:1);
			CS.psi_t(:,:,ks) = psi(:,:,ks) * Qt;
			
			% Normalize psi_t, s.t. integral(psi_t' * psi_t) = 1
			scfac = 1 ./ sqrt(sum(repmat(S.W,1,S.Ns_top) .* (CS.psi_t(:,:,ks) .* conj(CS.psi_t(:,:,ks))),1));
			CS.psi_t(:,:,ks) = bsxfun(@times, CS.psi_t(:,:,ks), scfac);
			
			% Normalize psi, s.t. integral(psi_new' * psi_new) = 1
			scfac = 1 ./ sqrt(sum(repmat(S.W,1,S.Nev) .* (psi(:,:,ks) .* conj(psi(:,:,ks))),1));
			% psi(:,:,ks) = psi(:,:,ks) * diag(scfac);
			psi(:,:,ks) = bsxfun(@times, psi(:,:,ks), scfac);
		else % use standard eigensolver to solve the subspace eigenproblem
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
				a0(ks) = min(S.EigVal(:,ks));
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
				a0(ks) = min(S.EigVal(:,ks));
			end
			% Set the filter cut off
			%S.lambda_f + log10(1e6-1)/S.bet + 0.5; 
			%lambda_cutoff(ks) = max(S.EigVal(:,ks)) + 0.10; 
		end
		
		if count == 1 && S.ForceCount == 1
			lambda_cutoff(ks) = 0.5 * (bup(ks) + a0(ks));
		else
			% Set the filter cut off
			%S.lambda_f + log10(1e6-1)/S.bet + 0.5; 
			lambda_cutoff(ks) = max(S.EigVal(:,ks)) + 0.10; 
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


