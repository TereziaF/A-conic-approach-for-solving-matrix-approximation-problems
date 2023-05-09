function[X,t_var,S,V] = OPP_CI1_bicri_L1(U,C,A,B,W,m,n,p,q,alpha)

cvx_begin
variable X(m,n) 
variable t_var(1,1)
variable S(p,q)
variable V(m+n, m+n) symmetric
minimize t_var + alpha*sum(diag(U*V))

V == [eye(m), X; X', eye(n)];

V == semidefinite(m+n);

-S <= W.*(C-A*X*B) <= S;

S'*ones(p,1) <= t_var*ones(q,1);

cvx_end
end