function[U] = CI2(V,m,n,k)

cvx_begin 
variable U(n+m,n+m) symmetric
minimize sum(diag(U*V))

sum(diag(U)) == n+m-k;
eye(n+m) - U == semidefinite(m+n);
U == semidefinite(n+m);

cvx_end