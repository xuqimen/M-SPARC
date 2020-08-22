function [Qt, Dt] = eigsolver_ChebFilt(Hs, Qt, m, n, a, b, a0)
% @brief Use Chebyshev filter to solve for the top Nt eigenvalues and eigenvectors
% of the matrix Hs.
% @param a0  Lower bound of lambda(-Hs), doesn't need to be very accurate.
% @param b   Upper bound of lambda(-Hs), needs to be >= actual max eig(-Hs).
% @param a   Cutoff of lambda(-Hs), near -1 * the Nt-th eigenvalue of Hs (with a
% little buffer to the right).

fprintf(2,'Filter bounds for eigsolver_ChebFilt: a = %f, b = %f, a0 = %f\n',a,b,a0);

% [Qt, Dt] = eigsolver_ChebFilt_test2(Hs, Qt, m, a, b, a0);
% return;

% [Qt, Dt] = eigsolver_ChebFilt_worked(Hs, Qt, m, n, a, b, a0);
% return;


[Qt, Dt] = eigsolver_ChebFilt_test3(Hs, Qt, m, n, a, b, a0);
return;

% n = 1;
for r = 1:n
	%Qt = chebyshev_filter_dense(-Hs,Qt,m,a,b,a0);
	Qt = chebyshev_filter_dense(-Hs,Qt,m,a,b,a0);
	Qt = orthChol(Qt);
end

Ht_s = Qt' * Hs * Qt;
Ht_s = 0.5 * (Ht_s + Ht_s');
[Qt_s, Dt] = eig(Ht_s);
Qt = Qt * Qt_s;

end


function [Qt, Dt] = eigsolver_ChebFilt_worked(Hs, Qt, m, n, a, b, a0)
% @brief Use Chebyshev filter to solve for the top Nt eigenvalues and eigenvectors
% of the matrix Hs.
% @param a0  Lower bound of lambda(-Hs), doesn't need to be very accurate.
% @param b   Upper bound of lambda(-Hs), needs to be >= actual max eig(-Hs).
% @param a   Cutoff of lambda(-Hs), near -1 * the Nt-th eigenvalue of Hs (with a
% little buffer to the right).

%n = 1;
for r = 1:n
	Qt = chebyshev_filter_dense(-Hs,Qt,m,a,b,a0);
	Qt = orthChol(Qt);
end

Ht_s = Qt' * Hs * Qt;
Ht_s = 0.5 * (Ht_s + Ht_s');
[Qt_s, Dt] = eig(Ht_s);

Qt = Qt * Qt_s;

end



function [Qt, Dt] = eigsolver_ChebFilt_test3(Hs, Qt, m, n, a, b, a0)
% @brief Use Chebyshev filter to solve for the top Nt eigenvalues and eigenvectors
% of the matrix Hs.
% @param a0  Lower bound of lambda(-Hs), doesn't need to be very accurate.
% @param b   Upper bound of lambda(-Hs), needs to be >= actual max eig(-Hs).
% @param a   Cutoff of lambda(-Hs), near -1 * the Nt-th eigenvalue of Hs (with a
% little buffer to the right).

%n = 2;
for r = 1:n
	%Qt = chebyshev_filter_dense(-Hs,Qt,m,a,b,a0);
	Qt = chebyshev_filter_dense(Hs,Qt,m,a,b,a0);
	Qt = orthChol(Qt);
end

Ht_s = Qt' * Hs * Qt;
Ht_s = 0.5 * (Ht_s + Ht_s');
[Qt_s, Dt] = eig(Ht_s);
Qt = Qt * Qt_s;

	
% 	m
% 	Ht_s
% 	diag(Dt)
% 	Qt_s
% 	Qt
% 	xxx

end

