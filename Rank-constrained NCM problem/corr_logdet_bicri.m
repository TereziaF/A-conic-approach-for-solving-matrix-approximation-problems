function[X,hodnost,hodnost_final,g,g_final,norma,norm_final,cas,t,s,empirical_epsilon,Z,Y] = corr_logdet_bicri(X0,C,k,epsilon,M,alpha)

%inputs:
% C - empirical correlation matrix
% k - desired rank
% epsilon - tolerance
% M - maximum number of the same iterations

%dimension
n = size(C,1);

tic;
%inicialize
X = X0;
hodnost = sum(eig(X)>epsilon);
g = Inf;
norma = norm(C-X,'fro');
delta = 0.01;

%set counter of iterations
t=0;

%set counter of the same iterations
s=0;

while (hodnost(end) > k && s < M)

t = t+1;

inverzna_matica = inv(X + delta*eye(n));

cvx_begin quiet
variable X(n,n) symmetric
variable Z(n,n)
variable Y(2*n, 2*n) symmetric
minimize sum(diag(Z)) + alpha*sum(diag(inverzna_matica*X))

for i=1:n
    X(i,i) == 1;
end

Y == [eye(n), (C-X)'; C-X, Z];

Y == semidefinite(2*n);

X == semidefinite(n);

cvx_end

hodnost = [hodnost;sum(eig(X)>epsilon)];
g = [g; sum(diag(Z))];
norma = [norma; norm(C-X,'fro')];

if hodnost(end) == hodnost(end-1)
    s = s+1;
else
    s = 1;
end

end

vh = eig(X);
empirical_epsilon = vh(n-k);
g_final = sum(diag(Z));
norm_final = norm(C-X,'fro');
hodnost_final = hodnost(end);
cas = toc;

end