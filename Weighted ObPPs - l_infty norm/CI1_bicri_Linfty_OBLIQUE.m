function[X,t_var,S,V,G] = CI1_bicri_Linfty_OBLIQUE(U,C,A,B,W,m,n,p,q,alpha)

cvx_begin
variable X(m,n) 
variable G(n,n) symmetric
variable t_var(1,1)
variable S(p,q)
variable V(m+n, m+n) symmetric
minimize t_var + alpha*sum(diag(U*V))

%INSERT ADDITIONAL LINEAR OR SEMIDEFINITE CONSTRAINTS HERE

diag(G) == ones(n,1);

V == [eye(m), X; X', G];

V == semidefinite(m+n);

-S <= W.*(C-A*X*B) <= S;

S*ones(q,1) <= t_var*ones(p,1);

cvx_end
end
