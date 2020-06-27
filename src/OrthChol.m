function Q = OrthChol(V)
% ORTHCHOL uses cholesky decomposition to orthogonalize V

R = chol(V'*V);
tic
Q = V*inv(R);
toc