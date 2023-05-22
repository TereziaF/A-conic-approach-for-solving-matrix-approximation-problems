function[X,hodnost,hodnost_final,g,g_final,norma,norm_final,cas,t,s,empirical_epsilon,sum_eig_final,sum_eig_real,t_var,S,V] = OPP_cvx_iter_bicri_Linfty(C,A,B,W,k,epsilon,M,alpha)

%inputs:
% C, A, B - data of the problem
% W - matrix specifying missing elements of C
% k - desired rank
% epsilon - tolerance
% M - maximum number of the same iterations
% alpha - relative weight for bi-criterion problem

%outputs:
% X - orthogonal solution
% hodnost - rank of variable V in all iterations
% hodnost_final - rank of solution V
% g - values of the objective in all iterations
% g_final - optimal value of the reformulated problem
% norma - values of the original objective in all iterations
% norm_final - optimal value of the original problem
% cas - computation time
% t - number of iterations
% s - number of "same-rank" iterations
% empirical_epsilon - value of the k-largest eigenvalue
% sum_eig_final - optimal value of OPP_CI2 in final iteration
% sum_eig_real - sum of (m+n-k) smallest eigenvalues of solution V
% t_var,S,V - solutions of the reformulated problem

%dimensions
p = size(C,1);
q = size(C,2);
m = size(A,2);
n = size(B,1);
k=m;

tic;

%inicialize
[X,t_var,S,V] = SDPrelax_Linfty(A, B, C, W, m, n, p, q)
U = OPP_CI2(V,m,n,k);

hodnost = sum(eig(V)>epsilon);
sum_eig = sum(diag(V*U));
norma = norm(W.*(C-A*X*B),Inf);
g = t_var;

%set counter of iterations
t=0;

%set counter of the same iterations
s=0;

%algorithm
while (hodnost(end) > k && s < M)

t = t+1;

[X,t_var,S,V] = OPP_CI1_bicri_Linfty(U,C,A,B,W,m,n,p,q,alpha);
U = OPP_CI2(V,m,n,k);

%saving values
hodnost = [hodnost;sum(eig(V)>epsilon)];
sum_eig = [sum_eig;sum(diag(V*U))];
norma = [norma; norm(W.*(C-A*X*B),Inf)];
g = [g; t_var];

if hodnost(end) == hodnost(end-1)
    s = s+1;
else
    s = 1;
end

end

%specifying outputs
vh = eig(V);
empirical_epsilon = vh(n+m-k);
norm_final = norm(W.*(C-A*X*B),Inf);
g_final = t_var;
hodnost_final = hodnost(end);
sum_eig_final = sum_eig(end);
sum_eig_real = sum(vh(1:(n+m-k)));
cas = toc;

end
