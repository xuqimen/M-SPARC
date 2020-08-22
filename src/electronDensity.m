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

if S.CSFlag == 1
	S = electronDensity_CS(S);
	%S = electronDensity_CS_toprot(S); % use rotated top states and un-rotated bottom states
	%S = electronDensity_CS_unrot(S); % use un-rotated states
	return;
end

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
	
end


function S = electronDensity_CS(S)
% @brief    Calculate new density based on the new states using Complementary Subspace method.
S.rho = 0 * S.rho;
if S.nspin == 1
	for kpt =1:S.tnkpt
		S.rho = S.rho + 2*S.wkpt(kpt)* sum( bsxfun(@times,S.psi(:,:,kpt).*conj(S.psi(:,:,kpt)),ones(1,S.Nev)),2);
		S.rho = S.rho - 2*S.wkpt(kpt)* sum( bsxfun(@times,S.psi_t(:,:,kpt).*conj(S.psi_t(:,:,kpt)),1-S.occ(S.CS_index,kpt)'),2);
	end
	S.rho = real(S.rho);
else
	ks = 1;
	for spin =1:S.nspin
		for kpt =1:S.tnkpt
			S.rho(:,spin+1) = S.rho(:,spin+1) + S.wkpt(kpt)* sum( bsxfun(@times,S.psi(:,:,ks).*conj(S.psi(:,:,ks)),ones(1,S.Nev)),2);
			S.rho(:,spin+1) = S.rho(:,spin+1) - S.wkpt(kpt)* sum( bsxfun(@times,S.psi_t(:,:,ks).*conj(S.psi_t(:,:,ks)),1-S.occ(S.CS_index,ks)'),2);
			ks = ks + 1;
		end
	end
	S.rho(:,2:3) = real(S.rho(:,2:3));
	S.rho(:,1) = S.rho(:,2) + S.rho(:,3);
end

fprintf('int_rho = %.15f\n',sum(S.rho(:,1)) * S.dV);
	
end


function S = electronDensity_CS_unrot(S)
% @brief    Calculate new density based on the new states (unrotated states).

% use the unrotated psi for evaluating density, with the standard formula
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
	
end



function S = electronDensity_CS_toprot(S)
% @brief    Calculate new density based on the new states (unrotated states).

% use the unrotated bottom states and rotated top states for evaluating density, with the standard formula
S.psi(:,S.CS_index,:) = S.psi_t;

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
	
end