function[X_hat,hodnost_hat,g_hat,norma_hat,gamma,hodnost,iter,cas,hodnost_final_CI_MOD,g_final_CI_MOD,norm_final_CI_MOD,cas_CI_MOD,t_CI_MOD,s_CI_MOD,empirical_epsilon_CI_MOD,sum_eig_final,sum_eig_real] = bisection_CI_OBLIQUE_L2(g0,g1,X0,X1,V0,V1,C,A,B,W, k,epsilon,M)

%dimension
p = size(C,1);
q = size(C,2);
m = size(A,2);
n = size(B,1);
k=m;

a1 = tic;
%initialize
X_hat = X1;
V_hat = V1;
g_hat = 0;

%set counter of iterations
iter = 1;

hodnost_final_CI_MOD(iter) = sum(eig(V1)>epsilon);
g_final_CI_MOD(iter) = g1;
norm_final_CI_MOD(iter) = norm(W.*(C-A*X1*B),2);
cas_CI_MOD(iter) = 0;
t_CI_MOD(iter) = 0;
s_CI_MOD(iter) = 0;
vh = eig(V1);
empirical_epsilon_CI_MOD(iter) = vh(n+m-k);
hodnost(iter) = sum(eig(V1)>epsilon);
gamma(iter) = 0;
sum_eig_final(iter) = 0;
sum_eig_real(iter) = sum(vh(1:n+m-k));

%interval
lb = g0;
ub = g1;

while abs(ub-lb) > epsilon
    
    iter = iter + 1;
    
    gamma(iter) = (lb+ub)/2;
    
    [X,hodnost_CI_MOD,hodnost_final_CI_MOD(iter),g_CI_MOD,g_final_CI_MOD(iter),norma_CI_MOD,norm_final_CI_MOD(iter),cas_CI_MOD(iter),t_CI_MOD(iter),s_CI_MOD(iter),empirical_epsilon_CI_MOD(iter),sum_eig_final(iter),sum_eig_real(iter),s_var,Y,V,G] = cvx_iter_mod_OBLIQUE_L2(C,A,B,W,k,epsilon,M,gamma(iter));
    
    hodnost(iter) = sum(eig(V)>epsilon);
    
    if hodnost(iter) == k
        X_hat = X;
        V_hat = V;
        g_hat = s_var;
        ub = gamma(iter);
    else
        lb = gamma(iter);
    end
end

hodnost_hat = sum(eig(V_hat)>epsilon);
norma_hat = norm(W.*(C-A*X_hat*B),2);
a2 = toc;
cas = a2-a1;
end