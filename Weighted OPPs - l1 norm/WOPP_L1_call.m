%generate WOPP with L1 norm

%dimension
m = 10;
n = 3;
p = 20;
q = n;

%max. no. of iterations
maxiter = 10;
%tolerance constant
tol = 10^(-6);                                  

%generating solution
X_sol = RandOrthMat(m,n,1e-2);   

%randomly generated data
A = randi(30,p,m);               
B = eye(q);                     
%B = randi(30,n,q); 

%delta = 0;                      %zero gap
delta = randn(p,n);             %non-zero gap        
C = A*X_sol*B + 0.5*delta;                  
W = ones(p,q);

%SDP relaxation
tic;
[X_rel,t_rel,S_rel,V_rel] = SDPrelax_L1(A, B, C, W, m, n, p, q);
toc;

norm(C-A*X_rel, 1)                    %optimal value
norm(X_rel'*X_rel-eye(n),'fro')       %orthogonality criterion
sum(eig(V_rel)>tol)                   %rank


%logdet for RCCP with alpha > 0
alpha = 1;
[X_LOG_bicri,hodnost_LOG_bicri,hodnost_final_LOG_bicri,g_LOG_bicri,g_final_LOG_bicri,norma_LOG_bicri,norm_final_LOG_bicri,cas_LOG_bicri,t_LOG_bicri,s_LOG_bicri,empirical_epsilon_LOG_bicri,t_var_LOG_bicri,S_LOG_bicri,V_LOG_bicri] = OPP_logdet_bicri_L1(X_rel,V_rel,C,A,B,W,m,tol,maxiter,alpha);

norm_final_LOG_bicri                           %optimal value
hodnost_final_LOG_bicri                        %rank
norm(X_LOG_bicri'*X_LOG_bicri-eye(n),'fro')    %orthogonality criterion

%convex iteration for RCCP with alpha > 0
alpha = 100;
[X_CI_bicri,hodnost_CI_bicri,hodnost_final_CI_bicri,g_CI_bicri,g_final_CI_bicri,norma_CI_bicri,norm_final_CI_bicri,cas_CI_bicri,t_CI_bicri,s_CI_bicri,empirical_epsilon_CI_bicri,sum_eig_final_CI_bicri,sum_eig_real_CI_bicri,t_var_CI_bicri,S_CI_bicri,V_CI_bicri] = OPP_cvx_iter_bicri_L1(C,A,B,W,m,tol,maxiter,alpha);

norm_final_CI_bicri                            %optimal value
hodnost_final_CI_bicri                         %rank
norm(X_CI_bicri'*X_CI_bicri-eye(n),'fro')      %orthogonality criterion


%logdet for RCFP modified
gamma = t_rel;
[X_LOG_MOD,hodnost_LOG_MOD,hodnost_final_LOG_MOD,g_LOG_MOD,g_final_LOG_MOD,norma_LOG_MOD,norm_final_LOG_MOD,cas_LOG_MOD,t_LOG_MOD,s_LOG_MOD,empirical_epsilon_LOG_MOD,t_var_LOG_MOD,S_LOG_MOD] = OPP_logdet_mod_L1(X_rel,V_rel,C,A,B,W,m,tol,maxiter,gamma);

norm_final_LOG_MOD                             %optimal value
hodnost_final_LOG_MOD                          %rank
norm(X_LOG_MOD'*X_LOG_MOD-eye(n),'fro')        %orthogonality criterion

%convex iteration for RCFP modified
gamma = t_rel;
[X_CI_MOD,hodnost_CI_MOD,hodnost_final_CI_MOD,g_CI_MOD,g_final_CI_MOD,norma_CI_MOD,norm_final_CI_MOD,cas_CI_MOD,t_CI_MOD,s_CI_MOD,empirical_epsilon_CI_MOD,sum_eig_final_CI_MOD,sum_eig_real_CI_MOD,t_var_CI_MOD,S_CI_MOD] = OPP_cvx_iter_mod_L1(C,A,B,W,m,tol,maxiter,gamma);

norm_final_CI_MOD                             %optimal value
hodnost_final_CI_MOD                          %rank
norm(X_CI_MOD'*X_CI_MOD-eye(n),'fro')         %orthogonality criterion

%bisective algorithm with logdet for RCCP with alpha > 0
[X_hat,hodnost_hat,g_hat,norma_hat,gamma,hodnost,iter,cas,hodnost_final_LOG_MOD,g_final_LOG_MOD,norm_final_LOG_MOD,cas_LOG_MOD,t_LOG_MOD,s_LOG_MOD,empirical_epsilon_LOG_MOD] = OPP_bisection_logdet_L1(t_rel,g_final_LOG_bicri,X_rel,X_LOG_bicri,V_rel,V_LOG_bicri,C,A,B,W,m,tol,maxiter);

norma_hat                                %optimal value
hodnost_hat                              %rank
norm(X_hat'*X_hat-eye(n),'fro')          %orthogonality criterion

%bisective algorithm with convex iteration for RCCP with alpha > 0
[X_hat,hodnost_hat,g_hat,norma_hat,gamma,hodnost,iter,cas,hodnost_final_CI_MOD,g_final_CI_MOD,norm_final_CI_MOD,cas_CI_MOD,t_CI_MOD,s_CI_MOD,empirical_epsilon_CI_MOD] = OPP_bisection_CI_L1(t_rel,g_final_CI_bicri,X_rel,X_CI_bicri,V_rel,V_CI_bicri,C,A,B,W,m,tol,maxiter);

norma_hat                                %optimal value
hodnost_hat                              %rank
norm(X_hat'*X_hat-eye(n),'fro')          %orthogonality criterion