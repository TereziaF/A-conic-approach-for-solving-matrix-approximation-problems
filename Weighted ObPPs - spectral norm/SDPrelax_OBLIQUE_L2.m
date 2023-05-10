function[X,s,Y,V,G] = SDPrelax_OBLIQUE_L2(A, B, C, W, m, n, p, q)

cvx_begin
variable X(m,n) 
variable G(n,n) symmetric
variable s(1,1)
variable Y(p+q,p+q)
variable V(m+n, m+n) symmetric
minimize s

Y == [s*eye(p), W.*(C - A*X*B); W'.*(C - A*X*B)', s*eye(q)];

Y == semidefinite(q+p);

diag(G) == ones(n,1);

V == [eye(m), X; X', G];

V == semidefinite(m+n);

cvx_end

end