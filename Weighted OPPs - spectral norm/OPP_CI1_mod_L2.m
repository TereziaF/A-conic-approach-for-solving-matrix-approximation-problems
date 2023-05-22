function[X,Y,s,V] = OPP_CI1_mod_L2(U,C,A,B,W,m,n,p,q,gamma)

%inputs:
%U - output of OPP_CI2
%A,B,C - data of the given problem
%W - matrix specifying missing elements of C
%m,n,p,q - dimensions of data
%gamma - upper bound on the original objective

%outputs:
%X - orthogonal solution
%Y,s,V - solutions of the reformulated problem

cvx_begin
variable X(m,n) 
variable s(1,1)
variable Y(q+p, q+p) symmetric
variable V(m+n, m+n) symmetric
minimize sum(diag(U*V))

%INSERT ADDITIONAL LINEAR OR SEMIDEFINITE CONSTRAINTS HERE

Y == [s*eye(p), W.*(C - A*X*B); W'.*(C - A*X*B)', s*eye(q)];

Y == semidefinite(q+p);

V == [eye(m), X; X', eye(n)];

V == semidefinite(m+n);

 s <= gamma;

cvx_end
end
