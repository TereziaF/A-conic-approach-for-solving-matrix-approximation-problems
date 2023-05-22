function[X,Y,Z,V] = OPP_CI1_bicri(U,C,A,B,W,m,n,p,q,alpha)

%inputs:
%U - output of OPP_CI2
%A,B,C - data of the given problem
%W - matrix specifying missing elements of C
%m,n,p,q - dimensions of data
%alpha - relative weight for bi-criterion problem

%outputs:
%X - orthogonal solution
%Y,Z,V - solutions of the reformulated problem

cvx_begin sdp
variable X(m,n) 
variable Z(p,p) symmetric
variable Y(q+p, q+p) symmetric
variable V(m+n, m+n) symmetric
minimize sum(diag(Z)) + alpha*sum(diag(U*V))

%INSERT ADDITIONAL LINEAR OR SEMIDEFINITE CONSTRAINTS HERE

Y == [eye(q), W'.*(C - A*X*B)'; W.*(C - A*X*B), Z];

Y == semidefinite(q+p);

V == [eye(m), X; X', eye(n)];

V == semidefinite(m+n);

cvx_end
end