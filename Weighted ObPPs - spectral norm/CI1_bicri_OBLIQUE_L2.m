function[X,s_var,Y,V,G] = CI1_bicri_OBLIQUE_L2(U,C,A,B,W,m,n,p,q,alpha)

cvx_begin sdp
variable X(m,n) 
variable G(n,n) symmetric
variable Y(p+q,p+q) symmetric
variable V(m+n, m+n) symmetric
variable s_var(1,1)
minimize s_var + alpha*sum(diag(U*V))

Y == [s_var*eye(p), W.*(C - A*X*B); W'.*(C - A*X*B)', s_var*eye(q)];

Y == semidefinite(q+p);

diag(G) == ones(n,1);

V == [eye(m), X; X', G];

V == semidefinite(m+n);

cvx_end
end