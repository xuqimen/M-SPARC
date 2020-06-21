function Ci = ChebyshevCoeff(npl, AFUN, a, b)
% ChebyshevCoeff calculates the Chebyshev expansion coeffcients Ci s.t.
% AFUN(x) = sum_{i=0}^{npl} Ci * Ti(y), where y = (x-c)/e \in (-1,1), 
% c = (a+b)/2, e = (b-a)/2. Ti is the ith order Chebyshev polynomial of 
% the first kind.
% The standard Chebyshev expansion only works within (-1,1), since the 
% Chebyshev polynomial grow exponentially outside the (-1,1) range. This
% function tries to fit the funtion AFUN in the range of (a,b) by mapping
% the funtion from (a,b) to (-1,1). The coefficients calculated can be used
% to evaluate AFUN within (a,b) by 
%  AFUN(x) = sum_{i=0}^{npl} Ci * Ti(y), x \in (a,b), y = (x-c)/e \in (-1,1).

c = (b + a)/2;
e = (b - a)/2;
AFUN_map = @(y) AFUN(y*e+c);
Ci = ChebyshevCoeff_std(npl, AFUN_map);

end

function Ci = ChebyshevCoeff_std(npl, AFUN)
% ChebyshevCoeff calculates the Chebyshev expansion coeffcients Ci s.t.
% AFUN(x) = sum_{i=0}^{npl} Ci * Ti(x), where Ti is the ith order Chebyshev
% polynomial of the first kind.
% Warning: this function currently uses discrete orthogonality property of 
% Chebyshev polynomials which only give "approximate" coefficients. so if 
% the function is not smooth enough the coefficients will not be accurate 
% enough. so this function has to be replaced with its continuous counterpart 
% which evaluates oscillatory integrals and uses fft/fct.

d = zeros(1,npl+1);
Ci = zeros(1,npl+1);
for k = 0:npl
	y = cos(pi*(k-0.5+1)/(npl+1));
	d(k+1) = AFUN(y);
end

fac = 2.0/(npl+1);

for j = 0:npl
	sum = 0.0;
	for k = 0:npl
		sum = sum + d(k+1)*cos((pi*(j-1+1))*((k-0.5+1)/(npl+1)));
	end
	Ci(j+1) = fac*sum;
end

Ci(1) = Ci(1)/2.0;

end