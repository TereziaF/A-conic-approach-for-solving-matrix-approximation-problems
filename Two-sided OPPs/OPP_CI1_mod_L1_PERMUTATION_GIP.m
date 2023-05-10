function[X,t_var,S,V] = OPP_CI1_mod_L1_PERMUTATION_GIP(U,C,A,m,n,p,q,gamma)

cvx_begin
variable X(m,n) 
variable t_var(1,1)
variable S(p,q)
variable V(m+n, m+n) symmetric
minimize sum(diag(U*V))

V == [eye(m), X; X', eye(n)];

V == semidefinite(m+n);

-S <= X*C-A*X <= S;

S'*ones(p,1) <= t_var*ones(q,1);

t_var <= gamma;

for i=1:m
  for j=1:n
    X(i,j) >= 0;
  end
end


cvx_end
end