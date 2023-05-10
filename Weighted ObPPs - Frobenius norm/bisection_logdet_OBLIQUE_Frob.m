function[X_hat,hodnost_hat,g_hat,norma_hat,gamma,hodnost,iter,cas,hodnost_final_LOG_MOD,g_final_LOG_MOD,norm_final_LOG_MOD,cas_LOG_MOD,t_LOG_MOD,s_LOG_MOD,empirical_epsilon_LOG_MOD] = bisection_logdet_OBLIQUE_Frob(g0,g1,X0,X1,V0,V1,C,A,B,W,k,epsilon,M)


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
norm_final_LOG_MOD(iter) = norm(W.*(C-A*X1*B),'fro');
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

while abs(ub-lb) > epsilon
    
    iter = iter + 1;
    
    gamma(iter) = (lb+ub)/2;
    
    [X,hodnost_LOG_MOD,hodnost_final_LOG_MOD(iter),g_LOG_MOD,g_final_LOG_MOD(iter),norma_LOG_MOD,norm_final_LOG_MOD(iter),cas_LOG_MOD(iter),t_LOG_MOD(iter),s_LOG_MOD(iter),empirical_epsilon_LOG_MOD(iter),V,G,Y,Z] = logdet_mod_OBLIQUE_Frob(X0,V0,C,A,B,W,m,epsilon,M,gamma(iter));
    
    hodnost(iter) = sum(eig(V)>epsilon);
    
    if hodnost(iter) == k
        X_hat = X;
        V_hat = V;
        g_hat = sum(diag(Z));
        ub = gamma(iter);
    else
        lb = gamma(iter);
    end
end

hodnost_hat = sum(eig(V_hat)>epsilon);
norma_hat = norm(W.*(C-A*X_hat*B),'fro');
a2 = toc;
cas = a2-a1;
end