function[X,t_var,S,V] = OPP_CI1_bicri_L1(U,C,A,B,W,m,n,p,q,alpha)
  
%inputs:
%U - output of OPP_CI2
%A,B,C - data of the given problem
%W - matrix specifying missing elements of C
%m,n,p,q - dimensions of data
%alpha - relative weight for bi-criterion problem

%outputs:
%X - orthogonal solution
%Y,Z,V - solutions of the reformulated problem

cvx_begin
variable X(m,n) 
variable t_var(1,1)
variable S(p,q)
variable V(m+n, m+n) symmetric
minimize t_var + alpha*sum(diag(U*V))
  
%INSERT ADDITIONAL LINEAR OR SEMIDEFINITE CONSTRAINTS HERE

V == [eye(m), X; X', eye(n)];

V == semidefinite(m+n);

-S <= W.*(C-A*X*B) <= S;

S'*ones(p,1) <= t_var*ones(q,1);

cvx_end
end
