function[X,Y,s,V] = OPP_CI1_bicri_L2(U,C,A,B,W,m,n,p,q,alpha)

cvx_begin
variable X(m,n) 
variable s(1,1)
variable Y(q+p, q+p) symmetric
variable V(m+n, m+n) symmetric
minimize s + alpha*sum(diag(U*V))

Y == [s*eye(p), W.*(C - A*X*B); W'.*(C - A*X*B)', s*eye(q)];

Y == semidefinite(q+p);

V == [eye(m), X; X', eye(n)];

V == semidefinite(m+n);

cvx_end
end