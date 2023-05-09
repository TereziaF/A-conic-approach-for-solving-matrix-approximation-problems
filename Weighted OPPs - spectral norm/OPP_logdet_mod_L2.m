function[X,hodnost,hodnost_final,g,g_final,norma,norm_final,cas,t,s,empirical_epsilon,s_var,Y,V] = OPP_logdet_mod_L2(X0,V0,C,A,B,W,k,epsilon,M,gamma)

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
norma = norm(W.*(C-A*X*B),2);
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
variable s_var(1,1)
variable Y(q+p, q+p) symmetric
variable V(m+n, m+n) symmetric
minimize sum(diag(inverzna_matica*V))

Y == [s_var*eye(p), W.*(C - A*X*B); W'.*(C - A*X*B)', s_var*eye(q)];

Y == semidefinite(q+p);

V == [eye(m), X; X', eye(n)];

V == semidefinite(m+n);

s_var <= gamma;

cvx_end

hodnost = [hodnost;sum(eig(V)>epsilon)];
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
g_final = g(end);
norm_final = norm(W.*(C-A*X*B),2);
hodnost_final = hodnost(end);
cas = toc;

end