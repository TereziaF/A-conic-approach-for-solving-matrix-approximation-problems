function[X,hodnost,hodnost_final,g,g_final,norma,norm_final,cas,t,s,empirical_epsilon,t_var,S,V,G] = logdet_mod_Linfty_OBLIQUE(X0,V0,C,A,B,W,k,epsilon,M,gamma)

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
variable G(n,n) symmetric
variable t_var(1,1)
variable S(p,q)
variable V(m+n, m+n) symmetric
minimize sum(diag(inverzna_matica*V))

diag(G) == ones(n,1);

V == [eye(m), X; X', G];

V == semidefinite(m+n);

-S <= W.*(C-A*X*B) <= S;

S*ones(q,1) <= t_var*ones(p,1);

t_var <= gamma;

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
g_final = g(end);
norm_final = norm(W.*(C-A*X*B),Inf);
hodnost_final = hodnost(end);
cas = toc;

end