function[X_hat,hodnost_hat,g_hat,norma_hat,gamma,hodnost,iter,cas,hodnost_final_LOG_MOD,g_final_LOG_MOD,norm_final_LOG_MOD,cas_LOG_MOD,t_LOG_MOD,s_LOG_MOD,empirical_epsilon_LOG_MOD] = bisection_logdet(g0,g1,X0,X1,C,k,epsilon,M)

%inputs:
% C - empirical correlation matrix
% k - desired rank
% epsilon - tolerance
% M - maximum number of the same iterations

%dimension
n = size(C,1);

tic;
%initialize
X_hat = X1;

%set counter of iterations
iter = 1;

hodnost_final_LOG_MOD(iter) = sum(eig(X1)>epsilon);
g_final_LOG_MOD(iter) = g1;
norm_final_LOG_MOD(iter) = norm(C-X1,'fro');
cas_LOG_MOD(iter) = 0;
t_LOG_MOD(iter) = 0;
s_LOG_MOD(iter) = 0;
vh = eig(X1);
empirical_epsilon_LOG_MOD(iter) = vh(n-k);
hodnost(iter) = sum(eig(X1)>epsilon);
gamma(iter) = 0;

%interval
lb = g0;
ub = g1;

while abs(ub-lb) > epsilon
    
    iter = iter + 1;
    
    gamma(iter) = (lb+ub)/2;
    
    [X,hodnost_LOG_MOD,hodnost_final_LOG_MOD(iter),g_LOG_MOD,g_final_LOG_MOD(iter),norma_LOG_MOD,norm_final_LOG_MOD(iter),cas_LOG_MOD(iter),t_LOG_MOD(iter),s_LOG_MOD(iter),empirical_epsilon_LOG_MOD(iter),Z] = corr_logdet_mod(X0,C,k,epsilon,M,gamma(iter));
    
    hodnost(iter) = sum(eig(X)>epsilon);
    
    if hodnost(iter) == k
        X_hat = X;
        g_hat = sum(diag(Z));
        ub = gamma(iter);
    else
        lb = gamma(iter);
    end
end

hodnost_hat = sum(eig(X_hat)>epsilon);
norma_hat = norm(C-X_hat,'fro');
cas = toc;
end