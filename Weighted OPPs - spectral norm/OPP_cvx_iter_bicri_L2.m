function[X,hodnost,hodnost_final,g,g_final,norma,norm_final,cas,t,s,empirical_epsilon,sum_eig_final,sum_eig_real,s_var,Y,V] = OPP_cvx_iter_bicri_L2(C,A,B,W,k,epsilon,M,alpha)


%dimensions
p = size(C,1);
q = size(C,2);
m = size(A,2);
n = size(B,1);
k=m;

tic;

%inicialize
[X,Y,s_var,V] = SDPrelax_L2(A, B, C, W, m, n, p, q)
U = OPP_CI2(V,m,n,k);

hodnost = sum(eig(V)>epsilon);
sum_eig = sum(diag(V*U));
g = s_var;
norma = norm(W.*(C-A*X*B),2);

%set counter of iterations
t=0;

%set counter of the same iterations
s=0;

while (hodnost(end) > k && s < M)

t = t+1;

[X,Y,s_var,V] = OPP_CI1_bicri_L2(U,C,A,B,W,m,n,p,q,alpha);
U = OPP_CI2(V,m,n,k);

hodnost = [hodnost;sum(eig(V)>epsilon)];
sum_eig = [sum_eig;sum(diag(V*U))];
g = [g; s_var];
norma = [norma; norm(W.*(C-A*X*B),2)];

if hodnost(end) == hodnost(end-1)
    s = s+1;
else
    s = 1;
end

end

vh = eig(V);
empirical_epsilon = vh(n+m-k);
g_final = s_var;
norm_final = norm(W.*(C-A*X*B),2);
hodnost_final = hodnost(end);
sum_eig_final = sum_eig(end);
sum_eig_real = sum(vh(1:(n+m-k)));
cas = toc;

end