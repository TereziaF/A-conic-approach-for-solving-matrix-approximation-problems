function[X,hodnost,hodnost_final,g,g_final,norma,norm_final,cas,t,s,empirical_epsilon,s_var,Y,V] = OPP_logdet_bicri_L2(X0,V0,C,A,B,W,k,epsilon,M,alpha)

%inputs:
% X0 - solution X of SDP relaxation
% V0 - solution V of SDP relaxation
% C, A, B - data of the problem
% W - matrix specifying missing elements of C
% k - desired rank
% epsilon - tolerance
% M - maximum number of the same iterations
% alpha - relative weight for bi-criterion problem

%outputs:
% X - orthogonal solution
% hodnost - rank of variable V in all iterations
% hodnost_final - rank of solution V
% g - values of the objective in all iterations
% g_final - optimal value of the reformulated problem
% norma - values of the original objective in all iterations
% norm_final - optimal value of the original problem
% cas - computation time
% t - number of iterations
% s - number of "same-rank" iterations
% empirical_epsilon - value of the k-largest eigenvalue
% s_var,Y,V - solutions of the reformulated problem
    
%dimension
p = size(C,1);
q = size(C,2);
m = size(A,2);
n = size(B,1);

tic;
%inicialize
X = X0;
V = V0;
hodnost = sum(eig(V)>epsilon);
g = Inf;
norma = norm(W.*(C-A*X*B),2);
delta = 0.01;

%set counter of iterations
t=0;

%set counter of the same iterations
s=0;

%algorithm
while (hodnost(end) > k && s < M)

t = t+1;

inverzna_matica = inv(V + delta*eye(m+n));

cvx_begin quiet
variable X(m,n) 
variable s_var(1,1)
variable Y(q+p, q+p) symmetric
variable V(m+n, m+n) symmetric
minimize s_var + alpha*sum(diag(inverzna_matica*V))

%INSERT ADDITIONAL LINEAR OR SEMIDEFINITE CONSTRAINTS HERE

Y == [s_var*eye(p), W.*(C - A*X*B); W'.*(C - A*X*B)', s_var*eye(q)];

Y == semidefinite(q+p);

V == [eye(m), X; X', eye(n)];

V == semidefinite(m+n);

cvx_end

%saving values
hodnost = [hodnost;sum(eig(V)>epsilon)];
g = [g; s_var];
norma = [norma; norm(W.*(C-A*X*B),2)];

if hodnost(end) == hodnost(end-1)
    s = s+1;
else
    s = 1;
end

end

%specifying outputs
vh = eig(V);
empirical_epsilon = vh(n+m-k);
g_final = s_var;
norm_final = norm(W.*(C-A*X*B),2);
hodnost_final = hodnost(end);
cas = toc;

end
