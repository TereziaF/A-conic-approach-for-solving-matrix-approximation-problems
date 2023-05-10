function[X,V,G,Y,Z] = CI1_bicri_OBLIQUE_Frob(U,C,A,B,W,m,n,p,q,alpha)

cvx_begin sdp
variable X(m,n) 
variable G(n,n) symmetric
variable Y(p+q,p+q) symmetric
variable V(m+n, m+n) symmetric
variable Z(p,p) symmetric
minimize sum(diag(Z)) + alpha*sum(diag(U*V))

Y == [eye(q), W'.*(C - A*X*B)'; W.*(C - A*X*B), Z];

Y == semidefinite(q+p);

diag(G) == ones(n,1);

V == [eye(m), X; X', G];

V == semidefinite(m+n);

cvx_end
end