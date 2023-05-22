function[U] = CI2(V,m,n,k)

%inputs:
%V - output of CI1_...
%m,n - dimensions of data
%k - desired rank

%outputs:
%U - direction matrix for next iteration of the algorithm

cvx_begin 
variable U(n+m,n+m) symmetric
minimize sum(diag(U*V))

sum(diag(U)) == n+m-k;
eye(n+m) - U == semidefinite(m+n);
U == semidefinite(n+m);

cvx_end
