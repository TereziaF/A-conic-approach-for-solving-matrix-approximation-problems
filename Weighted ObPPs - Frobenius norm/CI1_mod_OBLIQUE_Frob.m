function[X,V,G,Y,Z] = CI1_mod_OBLIQUE_Frob(U,C,A,B,W,m,n,p,q,gamma)

%inputs:
%U - output of CI2_...
%A,B,C - data of the given problem
%W - matrix specifying missing elements of C
%m,n,p,q - dimensions of data
%gamma - upper bound on the original objective

%outputs:
%X - orthogonal solution
%Y,Z,V,G - solutions of the reformulated problem

cvx_begin sdp
variable X(m,n) 
variable G(n,n) symmetric
variable Y(p+q,p+q) symmetric
variable V(m+n, m+n) symmetric
variable Z(p,p) symmetric
minimize sum(diag(U*V))

%INSERT ADDITIONAL LINEAR OR SEMIDEFINITE CONSTRAINTS HERE

Y == [eye(q), W'.*(C - A*X*B)'; W.*(C - A*X*B), Z];

Y == semidefinite(q+p);

diag(G) == ones(n,1);

V == [eye(m), X; X', G];

V == semidefinite(m+n);

 sum(diag(Z)) <= gamma;

cvx_end
end
