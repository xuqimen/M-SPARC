function S = occupations_sq(S,ChebComp)
	
FermiEnergyEvaluator = @(lambda_f_g) fermiCalc(lambda_f_g,S,ChebComp);

S.lambda_f = fzero(FermiEnergyEvaluator,0);

%fprintf(2,' ------------------\n');
fprintf(' Fermi energy = %f\n',S.lambda_f);

% Calculate occupations
% if S.elec_T_type == 0 % fermi-dirac smearing
% 	S.occ = 1./(1+exp(S.bet*(S.EigVal-S.lambda_f)));
% elseif S.elec_T_type == 1 % gaussian smearing
% 	S.occ = 0.5 * (1.0 - erf(S.bet*(S.EigVal-S.lambda_f)));
% end

end
	
function f = fermiCalc(lambda_f_g, S, ChebComp)
	sq_npl = S.sq_npl;
	tr_Ti_vec = [ChebComp(:).tr_Ti];
	if S.elec_T_type == 0 % fermi-dirac smearing
		AFUN = @(x) 1./(1+exp(S.bet*(x-lambda_f_g)));
	else
		AFUN = @(x) 0.5*(1.0-erf(S.bet*(x-lambda_f_g)));
	end
	f = 0;
	ks = 1;
	for spin = 1:S.nspin
		for kpt = 1:S.tnkpt
			% for a given lambda_f_g, recalculate the coefficients cj
			c = ChebyshevCoeff(sq_npl, AFUN, ChebComp(1).eigmin(ks), ChebComp(1).eigmax(ks));
			if S.direct_densmat == 0
				f = f + S.occfac * S.wkpt(kpt) * dot(c,tr_Ti_vec(:));
			else
				fprintf(2,'USING DIRECT DENSITY MAT EVALUATION\n');
				f = f + S.occfac * S.wkpt(kpt) * trace(AFUN_MAT(AFUN,ChebComp(1).Hs));
			end
			ks = ks + 1;
		end
	end
	%f = f - S.Nelectron;
	f = f + S.NegCharge;
end

function fA = AFUN_MAT(fun,A)
	[C, D] = eig(A);
	fA = C * diag(fun(diag(D))) * C';
end
