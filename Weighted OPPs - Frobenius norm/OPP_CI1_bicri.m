function[X,Y,Z,V] = OPP_CI1_bicri(U,C,A,B,W,m,n,p,q,alpha)

cvx_begin sdp
variable X(m,n) 
variable Z(p,p) symmetric
variable Y(q+p, q+p) symmetric
variable V(m+n, m+n) symmetric
minimize sum(diag(Z)) + alpha*sum(diag(U*V))

Y == [eye(q), W'.*(C - A*X*B)'; W.*(C - A*X*B), Z];

Y == semidefinite(q+p);

V == [eye(m), X; X', eye(n)];

V == semidefinite(m+n);

cvx_end
end