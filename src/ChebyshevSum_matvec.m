function Y = ChebyshevSum_matvec(c,A,X,a,b)
% ChebyshevSum_matvec evaluates Y = p(A)*X, where p is an expansion of the
% Chebyshev polynomials p(x) = sum_{i=0}^n c(i) * Ti(x,a,b), where Ti(x,a,b)
% is Ti(x) with x in [a,b] mapped to [-1,1].

% this is a brute force way to evaluate the Chebyshev sum O(N^2)
% Y = c(1) * X;
% for i = 2:length(c)
% 	Y = Y + c(i) * Chebyshev_matvec(A,X,i-1,a,b);
% end

Y = ChebyshevSum_matvec_OnTheFly(c,A,X,a,b); % O(N)

% this is another implementation of the chebyshev sum using Clenshaw's
% method, check out the modified version of Clenshaw's method too, when x =
% +/-1 or close to +/-1.
% Input: x; c0, c1, ..., cN.
% Output:SN(x) = sum_{k=0}^N ck * Tk(x).
% bN+1 = 0; bN = cN.
% DO r = N - 1,N - 2, ..., 1:
%	br = 2 x br+1 - br+2 + cr.
% SN(x) = xb1 - b2 + c0.

end

function Ysum = ChebyshevSum_matvec_OnTheFly(C,A,X,a,b)
% ChebyshevSum_matvec evaluates Y = p(A)*X, where p is an expansion of the
% Chebyshev polynomials p(x) = sum_{i=0}^npl c(i) * Ti(x,a,b), where Ti(x,a,b)
% is Ti(x) with x in [a,b] mapped to [-1,1].

npl = length(C) - 1;
a0 = a;
e = (b-a)/2;
c = (b+a)/2;
%sigma = e/(a0 - c); 
sigma = e/(c - a0); 
sigma1 = sigma;
gamma = 2/sigma1;

Ysum = C(1) * X;
if (npl <= 0), return; end
%A(1:size(A,1)+1:end) = A(1:size(A,1)+1:end) - c;
%Y = (sigma1/e)*(A*X); % T1(A,a,b) * X
Y = (sigma1/e)*(A*X-c*X); % T1(A,a,b) * X
Ysum = Ysum + C(2)*Y;
ii = 2;
while(ii <= npl)
	sigma2 = 1/(gamma - sigma);
	% AX = A * Y;
	Ynew = (2*sigma2/e)*(A*Y - c*Y) - ((sigma*sigma2)*X);
	%Ynew = (2*sigma2/e)*(A*Y) - ((sigma*sigma2)*X);
	Ysum = Ysum + C(ii+1) * Ynew;
	X = Y;
	Y = Ynew;
	sigma = sigma2;
	ii = ii + 1;
end

end



% function Y = Chebyshev_matvec(A,X,m,a,b)
% % Chebyshev_matvec evaluates Y = Tm(A)*X, where Tm is the mth order 
% % Chebyshev polynomials of the first kind.
% Y = A*X;
% ii = 2;
% while(ii <= m)
% 	Ynew = 2*(A*Y) - X;
% 	X = Y;
% 	Y = Ynew;
% 	ii = ii + 1;
% end
% 
% end



function Y = Chebyshev_matvec(A,X,m,a,b)
% WARNING: this function acutally does Y = Tm(-A)*X!!! one way to fix is
% to switch a and b.
% Chebyshev_matvec evaluates Y = Tm(A)*X, where Tm is the mth order 
% Chebyshev polynomials of the first kind. The chebyshev polynomial is 
% scaled from the original [-1,1] interval to [a,b].
a0 = a;
e = (b-a)/2;
c = (b+a)/2;
%sigma = e/(a0 - c); 
sigma = e/(c - a0); 
sigma1 = sigma;
gamma = 2/sigma1;
%AX = A * X;
Y = (sigma1/e)*(A*X-c*X);

ii = 2;
while(ii <= m)
	sigma2 = 1/(gamma - sigma);
	% AX = A * Y;
	Ynew = (2*sigma2/e)*(A*Y - c*Y) - ((sigma*sigma2)*X);
	X = Y;
	Y = Ynew;
	sigma = sigma2;
	ii = ii + 1;
end

end
