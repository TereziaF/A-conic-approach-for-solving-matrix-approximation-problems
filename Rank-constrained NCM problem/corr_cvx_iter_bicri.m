function[X,hodnost,hodnost_final,g,g_final,norma,norm_final,cas,t,s,empirical_epsilon,sum_eig_final,sum_eig_real,Z,Y] = corr_cvx_iter_bicri(C,k,epsilon,M,alpha)

%inputs:
% C - empirical correlation matrix
% k - desired rank
% epsilon - tolerance
% M - maximum number of the same iterations

%dimension
n = size(C,1);

tic;

%inicialize
[X,bla1,bla2,bla3,bla4,bla5,Z] = corr_relaxation(C,k,epsilon);
U = corr_CI2(X,n,k);

hodnost = sum(eig(X)>epsilon);
sum_eig = sum(diag(X*U));
g = sum(diag(Z));
norma = norm(C-X,'fro');

%set counter of iterations
t=0;

%set counter of the same iterations
s=0;

while (hodnost(end) > k && s < M)

t = t+1;

[X,Y,Z] = corr_CI1_bicri(U,C,n,alpha);
U = corr_CI2(X,n,k);

hodnost = [hodnost;sum(eig(X)>epsilon)];
sum_eig = [sum_eig;sum(diag(X*U))];
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
sum_eig_final = sum_eig(end);
sum_eig_real = sum(vh(1:(n-k)));
cas = toc;

end