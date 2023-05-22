function[X,V,G,Y,Z] = SDPrelax_OBLIQUE_Frob(A, B, C, W, m, n, p, q)
 
%inputs:
%A,B,C - data of the given problem
%m,n,p,q - dimensions of data

%outputs:
%X - oblique solution
%Y,Z,V - solutions of the reformulated problem

cvx_begin
variable X(m,n) 
variable G(n,n) symmetric
variable Y(p+q,p+q) symmetric
variable V(m+n, m+n) symmetric
variable Z(p,p) symmetric
minimize sum(diag(Z))
 
%INSERT ADDITIONAL LINEAR OR SEMIDEFINITE CONSTRAINTS HERE

Y == [eye(q), W'.*(C - A*X*B)'; W.*(C - A*X*B), Z];

Y == semidefinite(q+p);

diag(G) == ones(n,1);

V == [eye(m), X; X', G];

V == semidefinite(m+n);

cvx_end

end
