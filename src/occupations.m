function S = occupations(S)

if S.CSFlag == 0 % standard method
	FermiEnergyEvaluator = @(lambda_f_g) fermiCalc(lambda_f_g,S);
else % Complementary Subspace
	FermiEnergyEvaluator = @(lambda_f_g) fermiCalc_CS(lambda_f_g,S);
end

S.lambda_f = fzero(FermiEnergyEvaluator,0);

%fprintf(2,' ------------------\n');
fprintf(' Fermi energy = %f\n',S.lambda_f);

% Calculate occupations
if S.CSFlag == 0 % standard method
	if S.elec_T_type == 0 % fermi-dirac smearing
		S.occ = 1./(1+exp(S.bet*(S.EigVal-S.lambda_f)));
	elseif S.elec_T_type == 1 % gaussian smearing
		S.occ = 0.5 * (1.0 - erf(S.bet*(S.EigVal-S.lambda_f)));
	end
else
	S.occ = ones(size(S.EigVal));
	if S.elec_T_type == 0 % fermi-dirac smearing
		S.occ(S.CS_index,:) = 1./(1+exp(S.bet*(S.EigVal(S.CS_index,:)-S.lambda_f)));
	elseif S.elec_T_type == 1 % gaussian smearing
		S.occ(S.CS_index,:) = 0.5 * (1.0 - erf(S.bet*(S.EigVal(S.CS_index,:)-S.lambda_f)));
	end
end
S.occ
%xxxx

end


	
function f = fermiCalc(lambda_f_g, S)
	f = 0;
	ks = 1;
	for spin = 1:S.nspin
		for kpt = 1:S.tnkpt
			if S.elec_T_type == 0 % fermi-dirac smearing
				f = f + S.occfac * sum(S.wkpt(kpt)./(1+exp(S.bet*(S.EigVal(:,ks)-lambda_f_g)))) ;
			else
				f = f + S.occfac * sum(S.wkpt(kpt)*0.5*(1.0-erf(S.bet*(S.EigVal(:,ks)-lambda_f_g)))) ;
			end
			ks = ks + 1;
		end
	end
	%f = f - S.Nelectron;
	f = f + S.NegCharge;
end


function f = fermiCalc_CS(lambda_f_g, S)
	f = 0;
	ks = 1;
	for spin = 1:S.nspin
		for kpt = 1:S.tnkpt
			if S.elec_T_type == 0 % fermi-dirac smearing
				f = f + S.occfac * sum(S.wkpt(kpt)./(1+exp(S.bet*(S.EigVal(S.CS_index,ks)-lambda_f_g)))) ;
			else
				f = f + S.occfac * sum(S.wkpt(kpt)*0.5*(1.0-erf(S.bet*(S.EigVal(S.CS_index,ks)-lambda_f_g)))) ;
			end
			ks = ks + 1;
		end
	end
	%f = f - S.Nelectron;
	f = f + S.NegCharge + S.occfac * (S.Nev - S.Ns_top);
end

