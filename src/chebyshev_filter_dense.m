function Y = chebyshev_filter_dense(H,X,m,a,b,a0)
% @brief    CHEBYSHEV_FILTER_DENSE performs chebyshev fitering on the given
%           states. This function takes a given dense matrix instead of
%           sparse mat-vec routine.
%
% @param X  Current states to be filtered.
% @param m  Chebyshev polynomial degree.
% @param a  Lower bound of the filter.
% @param b  Upper bound of the filter.
% @param a0 Scaling constant of the filter.
%
% @ref:  
% Y. Zhou, et al, 2016. Parallel Self-Consistent-Field Calculations via
% Chebyshev-Filtered Subspace Acceleration.
% https://www-users.cs.umn.edu/~saad/PDF/umsi-2006-101.pdf
%
% @authors  Qimen Xu <qimenxu@gatech.edu>
%           Abhiraj Sharma <asharma424@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech
%a0 = a; % TODO: remove after check
e = (b-a)/2;
c = (b+a)/2;
sigma = e/(a0 - c); sigma1 = sigma;
%sigma = -e/(a0 - c); sigma1 = sigma; % TODO: recover after test
gamma = 2/sigma1;
HX = H * X;
Y = (sigma1/e)*(HX-c*X);

ii = 2;
while(ii <= m)
	sigma2 = 1/(gamma - sigma);
	HX = H * Y;
	Ynew = (2*sigma2/e)*(HX - c*Y) - ((sigma*sigma2)*X);
	X = Y;
	Y = Ynew;
	sigma = sigma2;
	ii = ii + 1;
end

end
