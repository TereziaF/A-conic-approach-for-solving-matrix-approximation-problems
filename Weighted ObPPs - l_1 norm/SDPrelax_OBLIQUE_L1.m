function[X,t,S,V,G] = SDPrelax_OBLIQUE_L1(A, B, C, W, m, n, p, q)

cvx_begin
variable X(m,n) 
variable G(n,n) symmetric
variable t(1,1)
variable S(p,q)
variable V(m+n, m+n) symmetric
minimize t

%INSERT ADDITIONAL LINEAR OR SEMIDEFINITE CONSTRAINTS HERE

diag(G) == ones(n,1);

V == [eye(m), X; X', G];

V == semidefinite(m+n);

-S <= W.*(C-A*X*B) <= S;

S'*ones(p,1) <= t*ones(q,1);

cvx_end

end
