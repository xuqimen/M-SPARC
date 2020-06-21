function S = electronDensity(S)
% @brief    Calculate new density based on the new states.
%
% @param S  A struct that contains the relevant fields.
%
% @authors  Qimen Xu <qimenxu@gatech.edu>
%           Abhiraj Sharma <asharma424@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech
%

if S.densMatFlag == 0
	S.rho = 0 * S.rho;
	if S.nspin == 1
		for kpt =1:S.tnkpt
			S.rho = S.rho + 2*S.wkpt(kpt)* sum( bsxfun(@times,S.psi(:,:,kpt).*conj(S.psi(:,:,kpt)),S.occ(:,kpt)'),2);
		end
		S.rho = real(S.rho);
	else
		ks = 1;
		for spin =1:S.nspin
			for kpt =1:S.tnkpt
				S.rho(:,spin+1) = S.rho(:,spin+1) + S.wkpt(kpt)* sum( bsxfun(@times,S.psi(:,:,ks).*conj(S.psi(:,:,ks)),S.occ(:,ks)'),2);
				ks = ks + 1;
			end
		end
		S.rho(:,2:3) = real(S.rho(:,2:3));
		S.rho(:,1) = S.rho(:,2) + S.rho(:,3);
	end
else
	% use density matrix to evaluate rho
	S = electronDensity_sq(S, S.ChebComp);
end
	
end



function S = electronDensity_sq(S, ChebComp)
% @brief    Calculate new density based on the new density matrix.

% first evaluate the subspace density matrix using the correct fermi-level
S.Ds = SubDensMat(S, S.lambda_f, ChebComp);

% the density matrix is Psi * Ds * Psi', we just need the diagonal
% note here Psi has to be scaled based on L2 norm, but our S.Psi has l2
% norm, so we need to divide the diagonal by dV
S.rho = 0 * S.rho;
if S.nspin == 1
	for kpt =1:S.tnkpt
		S.rho = S.rho + 2*S.wkpt(kpt)* sum((S.psi(:,:,kpt)*S.Ds).*conj(S.psi(:,:,kpt)),2);
	end
	S.rho = real(S.rho);
else
	ks = 1;
	for spin =1:S.nspin
		for kpt =1:S.tnkpt
			S.rho(:,spin+1) = S.rho(:,spin+1) + S.wkpt(kpt)* sum((S.psi(:,:,ks)*S.Ds).*conj(S.psi(:,:,ks)),2);
			ks = ks + 1;
		end
	end
	S.rho(:,2:3) = real(S.rho(:,2:3));
	S.rho(:,1) = S.rho(:,2) + S.rho(:,3);
end

% scale the diagonal by dV
S.rho = S.rho / S.dV;

fprintf(2, 'int rho = %.15f\n', dot(S.rho(:,1),S.W));

% S.rho = 0 * S.rho;
% if S.nspin == 1
% 	for kpt =1:S.tnkpt
% 		S.rho = S.rho + 2*S.wkpt(kpt)* sum( bsxfun(@times,S.psi(:,:,kpt).*conj(S.psi(:,:,kpt)),S.occ(:,kpt)'),2);
% 	end
% 	S.rho = real(S.rho);
% else
% 	ks = 1;
% 	for spin =1:S.nspin
% 		for kpt =1:S.tnkpt
% 			S.rho(:,spin+1) = S.rho(:,spin+1) + S.wkpt(kpt)* sum( bsxfun(@times,S.psi(:,:,ks).*conj(S.psi(:,:,ks)),S.occ(:,ks)'),2);
% 			ks = ks + 1;
% 		end
% 	end
% 	S.rho(:,2:3) = real(S.rho(:,2:3));
% 	S.rho(:,1) = S.rho(:,2) + S.rho(:,3);
% end
	
end


function Ds = SubDensMat(S, lambda_f, ChebComp)

sq_npl = S.sq_npl;
if S.elec_T_type == 0 % fermi-dirac smearing
	AFUN = @(x) 1./(1+exp(S.bet*(x-lambda_f)));
else
	AFUN = @(x) 0.5*(1.0-erf(S.bet*(x-lambda_f)));
end

Ds = zeros(S.Nev,S.Nev);
ks = 1;
for spin = 1:S.nspin
	for kpt = 1:S.tnkpt
		a = ChebComp(1).eigmin(ks);
		b = ChebComp(1).eigmax(ks);
		c = ChebyshevCoeff(sq_npl, AFUN, a, b);
		% TODO: Tj(Hs) are already known! no need to calculate them again
		if S.direct_densmat == 0
			Ds = Ds + S.wkpt(kpt) * ChebyshevSum_matvec(c,ChebComp(1).Hs,eye(S.Nev),a,b); 
		else
			fprintf(2,'USING DIRECT DENSITY MAT EVALUATION\n');
			Ds = Ds + S.wkpt(kpt) *  AFUN_MAT(AFUN,ChebComp(1).Hs);
		end
	end
end

%Ds = (Ds + Ds') * 0.5; 
fprintf(2, ' trace of Ds is %.15f\n',trace(Ds));
%eig(Ds)
end

function fA = AFUN_MAT(fun,A)
	[C, D] = eig(A);
	fA = C * diag(fun(diag(D))) * C';
end