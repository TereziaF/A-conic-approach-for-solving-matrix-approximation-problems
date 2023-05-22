function[X,Y,s,V] = SDPrelax_L2(A, B, C, W, m, n, p, q)
  
%inputs:
%A,B,C - data of the given problem
%m,n,p,q - dimensions of data

%outputs:
%X - orthogonal solution
%Y,s,V - solutions of the reformulated problem

cvx_begin sdp
variable X(m,n) 
variable s(1,1)
variable Y(q+p, q+p) symmetric
variable V(m+n, m+n) symmetric
minimize s

%INSERT ADDITIONAL LINEAR OR SEMIDEFINITE CONSTRAINTS HERE

Y == [s*eye(p), W.*(C - A*X*B); W'.*(C - A*X*B)', s*eye(q)];

Y == semidefinite(q+p);

V == [eye(m), X; X', eye(n)];

V == semidefinite(m+n);

cvx_end

end
