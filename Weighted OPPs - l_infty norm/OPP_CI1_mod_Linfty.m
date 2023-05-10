function[X,t_var,S,V] = OPP_CI1_mod_Linfty(U,C,A,B,W,m,n,p,q,gamma)

cvx_begin
variable X(m,n) 
variable t_var(1,1)
variable S(p,q)
variable V(m+n, m+n) symmetric
minimize sum(diag(U*V))

V == [eye(m), X; X', eye(n)];

V == semidefinite(m+n);

-S <= W.*(C-A*X*B) <= S;

S*ones(q,1) <= t_var*ones(p,1);

t_var <= gamma;

cvx_end
end