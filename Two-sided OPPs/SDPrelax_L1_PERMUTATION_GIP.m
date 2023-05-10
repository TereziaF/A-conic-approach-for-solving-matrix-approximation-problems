function[X,t,S,V] = SDPrelax_L1_PERMUTATION_GIP(A,C,m,n,p,q)

cvx_begin
variable X(m,n) 
variable t(1,1)
variable S(p, q)
variable V(m+n, m+n) symmetric
minimize t

-S <= X*C-A*X <= S;

S'*ones(p,1) <= t*ones(q,1);

V == [eye(m), X; X', eye(n)];

V == semidefinite(m+n);


for i=1:m
  for j=1:n
    X(i,j) >= 0;
  end
end

cvx_end

end