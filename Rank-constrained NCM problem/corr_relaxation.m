function[X,hodnost_final,g_final,norm_final,cas,empirical_epsilon,Z,Y] = corr_relaxation(C,k,epsilon)

%inputs:
% C - empirical correlation matrix
% k - desired rank
% epsilon - tolerance constant

%dimension
n = size(C,1);

tic;
%relaxation
cvx_begin quiet
variable X(n,n) symmetric
variable Z(n,n)
variable Y(2*n, 2*n) symmetric
minimize sum(diag(Z))

for i=1:n
    X(i,i) == 1;
end

Y == [eye(n), (C-X)'; C-X, Z];

Y == semidefinite(2*n);

X == semidefinite(n);

cvx_end

vh = eig(X);
empirical_epsilon = vh(n-k);
g_final = sum(diag(Z));
norm_final = norm(C-X,'fro');
hodnost_final = sum(eig(X)>epsilon);
cas = toc;
end