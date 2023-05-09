function[X,hodnost,hodnost_final,g,g_final,norma,norm_final,cas,t,s,empirical_epsilon,t_var,S,V] = OPP_logdet_bicri_Linfty(X0,V0,C,A,B,W,k,epsilon,M,alpha)

%dimension
p = size(C,1);
q = size(C,2);
m = size(A,2);
n = size(B,1);

tic;
%inicialize
X = X0;
V = V0;
hodnost = sum(eig(V)>epsilon);
g = Inf;
norma = norm(W.*(C-A*X*B),Inf);
delta = 0.01;

%set counter of iterations
t=0;

%set counter of the same iterations
s=0;

while (hodnost(end) > k && s < M)

t = t+1;

inverzna_matica = inv(V + delta*eye(m+n));

cvx_begin quiet
variable X(m,n) 
variable t_var(1,1)
variable S(p, q)
variable V(m+n, m+n) symmetric
minimize t_var + alpha*sum(diag(inverzna_matica*V))

V == [eye(m), X; X', eye(n)];

V == semidefinite(m+n);

-S <= W.*(C-A*X*B) <= S;

S*ones(q,1) <= t_var*ones(p,1);

cvx_end

hodnost = [hodnost;sum(eig(V)>epsilon)];
norma = [norma; norm(W.*(C-A*X*B),Inf)];
g = [g; t_var];

if hodnost(end) == hodnost(end-1)
    s = s+1;
else
    s = 1;
end

end

vh = eig(V);
empirical_epsilon = vh(n+m-k);
norm_final = norm(W.*(C-A*X*B),Inf);
g_final = t_var;
hodnost_final = hodnost(end);
cas = toc;

end