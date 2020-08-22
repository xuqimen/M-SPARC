function Q = orthChol(V)
% ORTHCHOL uses cholesky decomposition to orthogonalize V
M = V' * V
eig(M)
R = chol(M)
tic
Q = V*inv(R);
toc




% function Q = orthChol(V)
% % ORTHCHOL uses cholesky decomposition to orthogonalize V
% 
% R = chol(V'*V);
% tic
% Q = V*inv(R);
% toc