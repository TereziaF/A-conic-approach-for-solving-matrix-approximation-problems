function[X,hodnost,hodnost_final,g,g_final,norma,norm_final,cas,t,s,empirical_epsilon,V,G,Y,Z] = logdet_mod_OBLIQUE_Frob(X0,V0,C,A,B,W,k,epsilon,M,gamma)

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
norma = norm(W.*(C-A*X*B),'fro');
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
variable Y(p+q,p+q) symmetric
variable V(m+n, m+n) symmetric
variable Z(p,p) symmetric
minimize sum(diag(inverzna_matica*V))

Y == [eye(q), W'.*(C - A*X*B)'; W.*(C - A*X*B), Z];

Y == semidefinite(q+p);

diag(G) == ones(n,1);

V == [eye(m), X; X', G];

V == semidefinite(m+n);

sum(diag(Z)) <= gamma;

cvx_end

hodnost = [hodnost;sum(eig(V)>epsilon)];
g = [g; sum(diag(Z))];
norma = [norma; norm(W.*(C-A*X*B),'fro')];

if hodnost(end) == hodnost(end-1)
    s = s+1;
else
    s = 1;
end

end

vh = eig(V);
empirical_epsilon = vh(n+m-k);
g_final = g(end);
norm_final = norm(W.*(C-A*X*B),'fro');
hodnost_final = hodnost(end);
cas = toc;

end