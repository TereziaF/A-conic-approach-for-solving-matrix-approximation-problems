function[X,hodnost,hodnost_final,g,g_final,norma,norm_final,cas,t,s,empirical_epsilon,t_var,S,V] = OPP_logdet_mod_PERMUTATION_GIP(X0,V0,C,A,k,epsilon,M,gamma)

%dimension
p = size(C,1);
q = size(C,2);
m = size(A,2);
n = size(A,1);

tic;
%inicialize
X = X0;
V = V0;
hodnost = sum(eig(full(V))>epsilon);
g = Inf;
norma = norm(X*C-A*X,1);
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
minimize sum(diag(inverzna_matica*V))

-S <= X*C-A*X <= S;

S'*ones(p,1) <= t_var*ones(q,1);

V == [eye(m), X; X', eye(n)];

V == semidefinite(m+n);

for i=1:m
  for j=1:n
    X(i,j) >= 0;
  end
end

t_var <= gamma;

cvx_end

hodnost = [hodnost;sum(eig(full(V))>epsilon)];
g = [g; t_var];
norma = [norma; norm(X*C - A*X,1)];

if hodnost(end) == hodnost(end-1)
  s = s+1;
else
  s = 1;
end

end

vh = eig(full(V));
empirical_epsilon = vh(n+m-k);
g_final = g(end);
norm_final = norm(X*C-A*X,1);
hodnost_final = hodnost(end);
cas = toc;

end