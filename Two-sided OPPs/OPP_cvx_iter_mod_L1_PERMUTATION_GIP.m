function[X,hodnost,hodnost_final,g,g_final,norma,norm_final,cas,t,s,empirical_epsilon,sum_eig_final,sum_eig_real,t_var,S,V] = OPP_cvx_iter_mod_L1_PERMUTATION_GIP(C,A,k,epsilon,M,gamma)

%dimensions
p = size(C,1);
q = size(C,2);
m = size(A,2);
n = size(A,1);
k=m;

tic;

%inicialize
[X,t_var,S,V] = SDPrelax_L1_PERMUTATION_GIP(A, C, m, n, p, q);
U = OPP_CI2(V,m,n,k);

hodnost = sum(eig(V)>epsilon);
sum_eig = sum(diag(V*U));
g = t_var;
norma = norm(X*C-A*X,1);

%set counter of iterations
t=0;

%set counter of the same iterations
s=0;

while (hodnost(end) > k && s < M)

t = t+1;

[X,t_var,S,V] = OPP_CI1_mod_L1_PERMUTATION_GIP(U,C,A,m,n,p,q,gamma);
U = OPP_CI2(V,m,n,k);

hodnost = [hodnost;sum(eig(V)>epsilon)];
sum_eig = [sum_eig;sum(diag(V*U))];
g = [g; t_var];
norma = [norma; norm(X*C-A*X,1)];

if hodnost(end) == hodnost(end-1)
    s = s+1;
else
    s = 1;
end

end

vh = eig(V);
empirical_epsilon = vh(n+m-k);
g_final = t_var;
norm_final = norm(X*C-A*X,1);
hodnost_final = hodnost(end);
sum_eig_final = sum_eig(end);
sum_eig_real = sum(vh(1:(n+m-k)));
cas = toc;

end