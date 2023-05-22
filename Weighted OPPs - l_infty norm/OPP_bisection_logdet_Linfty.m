function[X_hat,hodnost_hat,g_hat,norma_hat,gamma,hodnost,iter,cas,hodnost_final_LOG_MOD,g_final_LOG_MOD,norm_final_LOG_MOD,cas_LOG_MOD,t_LOG_MOD,s_LOG_MOD,empirical_epsilon_LOG_MOD] = OPP_bisection_logdet_Linfty(g0,g1,X0,X1,V0,V1,C,A,B,W,k,epsilon,M)

%inputs:
% g0 - optimal value of SDP relaxation
% g1 - optimal value of a heuristic
% X0 - solution X of SDP relaxation
% X1 - solution X of a heuristic
% V0 - solution V of SDP relaxation
% V1 - solution V of a heuristic
% C, A, B - data of the problem
% W - matrix specifying missing elements of C
% k - desired rank
% epsilon - tolerance
% M - maximum number of the same iterations

%outputs:
% X_hat - orthogonal solution X
% hodnost_hat - rank of solution V
% g_hat - optimal value of the reformulated problem
% norm_hat - optimal value of the original problem
% gamma - values of gamma in all iterations
% hodnost - rank of solutions in all iterations
% iter - number of iterations
% cas - computation time

%dimension
p = size(C,1);
q = size(C,2);
m = size(A,2);
n = size(B,1);
k=m;

a1=tic;
%initialize
X_hat = X1;
V_hat = V1;
g_hat = 0;

%set counter of iterations
iter = 1;

hodnost_final_LOG_MOD(iter) = sum(eig(V1)>epsilon);
g_final_LOG_MOD(iter) = g1;
norm_final_LOG_MOD(iter) = norm(W.*(C-A*X1*B),Inf);
cas_LOG_MOD(iter) = 0;
t_LOG_MOD(iter) = 0;
s_LOG_MOD(iter) = 0;
vh = eig(V1);
empirical_epsilon_LOG_MOD(iter) = vh(n+m-k);
hodnost(iter) = sum(eig(V1)>epsilon);
gamma(iter) = 0;

%interval
lb = g0;
ub = g1;

%algorithm
while abs(ub-lb) > epsilon
    
    iter = iter + 1;
    
    gamma(iter) = (lb+ub)/2;
    
    [X,hodnost_LOG_MOD,hodnost_final_LOG_MOD(iter),g_LOG_MOD,g_final_LOG_MOD(iter),norma_LOG_MOD,norm_final_LOG_MOD(iter),cas_LOG_MOD(iter),t_LOG_MOD(iter),s_LOG_MOD(iter),empirical_epsilon_LOG_MOD(iter),t_var,S,V] = OPP_logdet_mod_Linfty(X0,V0,C,A,B,W,m,epsilon,M,gamma(iter));
    
    hodnost(iter) = sum(eig(V)>epsilon);
    
    if hodnost(iter) == k
        X_hat = X;
        V_hat = V;
        g_hat = t_var;
        ub = gamma(iter);
    else
        lb = gamma(iter);
    end
end

%specifying outputs
hodnost_hat = sum(eig(V_hat)>epsilon);
norma_hat = norm(W.*(C-A*X_hat*B),Inf);
a2 = toc;
cas = a2-a1;
end
