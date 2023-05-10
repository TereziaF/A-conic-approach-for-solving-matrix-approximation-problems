function[X,hodnost,hodnost_final,g,g_final,norma,norm_final,cas,t,s,empirical_epsilon,sum_eig_final,sum_eig_real,t_var,S,V,G] = cvx_iter_bicri_L1_OBLIQUE(C,A,B,W,k,epsilon,M,alpha)

%inputs:
% C - empirical correlation matrix
% k - desired rank
% epsilon - tolerance
% M - maximum number of the same iterations

%dimensions
p = size(C,1);
q = size(C,2);
m = size(A,2);
n = size(B,1);
k=m;

tic;

%inicialize
[X,t_var,S,V] = SDPrelax_OBLIQUE_L1(A, B, C, W, m, n, p, q)
U = CI2(V,m,n,k);

hodnost = sum(eig(V)>epsilon);
sum_eig = sum(diag(V*U));
norma = norm(W.*(C-A*X*B),1);
g = t_var;

%set counter of iterations
t=0;

%set counter of the same iterations
s=0;

while (hodnost(end) > k && s < M)

t = t+1;

[X,t_var,S,V,G] = CI1_bicri_L1_OBLIQUE(U,C,A,B,W,m,n,p,q,alpha);
U = CI2(V,m,n,k);

hodnost = [hodnost;sum(eig(V)>epsilon)];
sum_eig = [sum_eig;sum(diag(V*U))];
norma = [norma; norm(W.*(C-A*X*B),1)];
g = [g; t_var];

if hodnost(end) == hodnost(end-1)
    s = s+1;
else
    s = 1;
end

end

vh = eig(V);
empirical_epsilon = vh(n+m-k);
norm_final = norm(W.*(C-A*X*B),1);
g_final = t_var;
hodnost_final = hodnost(end);
sum_eig_final = sum_eig(end);
sum_eig_real = sum(vh(1:(n+m-k)));
cas = toc;

end