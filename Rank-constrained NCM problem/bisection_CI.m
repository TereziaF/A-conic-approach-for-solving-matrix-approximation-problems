function[X_hat,hodnost_hat,g_hat,norma_hat,gamma,hodnost,iter,cas,hodnost_final_CI_MOD,g_final_CI_MOD,norm_final_CI_MOD,cas_CI_MOD,t_CI_MOD,s_CI_MOD,empirical_epsilon_CI_MOD,sum_eig_final,sum_eig_real] = bisection_CI(g0,g1,X0,X1,C,k,epsilon,M)

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

hodnost_final_CI_MOD(iter) = sum(eig(X1)>epsilon);
g_final_CI_MOD(iter) = g1;
norm_final_CI_MOD(iter) = norm(C-X1,'fro');
cas_CI_MOD(iter) = 0;
t_CI_MOD(iter) = 0;
s_CI_MOD(iter) = 0;
vh = eig(X1);
empirical_epsilon_CI_MOD(iter) = vh(n-k);
hodnost(iter) = sum(eig(X1)>epsilon);
gamma(iter) = 0;
sum_eig_final(iter) = 0;
sum_eig_real(iter) = sum(vh(1:n-k));

%interval
lb = g0;
ub = g1;

while abs(ub-lb) > epsilon
    
    iter = iter + 1;
    
    gamma(iter) = (lb+ub)/2;
    
    [X,hodnost_CI_MOD,hodnost_final_CI_MOD(iter),g_CI_MOD,g_final_CI_MOD(iter),norma_CI_MOD,norm_final_CI_MOD(iter),cas_CI_MOD(iter),t_CI_MOD(iter),s_CI_MOD(iter),empirical_epsilon_CI_MOD(iter),sum_eig_final(iter),sum_eig_real(iter),Z] = corr_cvx_iter_mod(C,k,epsilon,M,gamma(iter));
    
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