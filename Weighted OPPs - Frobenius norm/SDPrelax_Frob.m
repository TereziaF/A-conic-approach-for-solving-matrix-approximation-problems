function[X,Y,Z,V] = SDPrelax_Frob(A, B, C, W, m, n, p, q)

%inputs:
%A,B,C - data of the given problem
%m,n,p,q - dimensions of data

%outputs:
%X - orthogonal solution
%Y,Z,V - solutions of the reformulated problem

cvx_begin sdp
variable X(m,n) 
variable Z(p,p) symmetric
variable Y(q+p, q+p) symmetric
variable V(m+n, m+n) symmetric
minimize sum(diag(Z))

%INSERT ADDITIONAL LINEAR OR SEMIDEFINITE CONSTRAINTS HERE

Y == [eye(q), W'.*(C - A*X*B)'; W.*(C - A*X*B), Z];

Y == semidefinite(q+p);

V == [eye(m), X; X', eye(n)];

V == semidefinite(m+n);

cvx_end

end