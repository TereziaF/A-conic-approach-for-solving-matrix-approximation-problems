function[X,t,S,V] = SDPrelax_Linfty(A, B, C, W, m, n, p, q)

%inputs:
%A,B,C - data of the given problem
%m,n,p,q - dimensions of data

%outputs:
%X - orthogonal solution
%t,S,V - solutions of the reformulated problem

cvx_begin
variable X(m,n) 
variable t(1,1)
variable S(p, q)
variable V(m+n, m+n) symmetric
minimize t

%INSERT ADDITIONAL LINEAR OR SEMIDEFINITE CONSTRAINTS HERE

V == [eye(m), X; X', eye(n)];

V == semidefinite(m+n);

-S <= W.*(C-A*X*B) <= S;

S*ones(q,1) <= t*ones(p,1);

cvx_end

end
