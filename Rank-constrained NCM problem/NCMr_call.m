%generate empirical correlation matrix
n=5;
C = approx_corr(n);

%define inputs for methods
epsilon = 10^(-6);
k = 2;
alpha = 100;
M = 20;

%relaxation for RCCP
[X_REL,hodnost_final_REL,g_final_REL,norm_final_REL,cas_REL,empirical_epsilon_REL,Z_REL,Y_REL] = corr_relaxation(C,k,epsilon);

gamma = g_final_REL;

%logdet for RCFP
[X_LOG,hodnost_LOG,hodnost_final_LOG,g_LOG,g_final_LOG,norma_LOG,norm_final_LOG,cas_LOG,t_LOG,s_LOG,empirical_epsilon_LOG,Z_LOG,Y_LOG] = corr_logdet(X_REL,C,k,epsilon,M);

%convex iteration for RCFP
[X_CI,hodnost_CI,hodnost_final_CI,g_CI,g_final_CI,norma_CI,norm_final_CI,cas_CI,t_CI,s_CI,empirical_epsilon_CI,sum_eig_final_CI,sum_eig_real_CI,Z_CI,Y_CI] = corr_cvx_iter(C,k,epsilon,M);

%logdet for RCCP with alpha > 0
[X_LOG_bicri,hodnost_LOG_bicri,hodnost_final_LOG_bicri,g_LOG_bicri,g_final_LOG_bicri,norma_LOG_bicri,norm_final_LOG_bicri,cas_LOG_bicri,t_LOG_bicri,s_LOG_bicri,empirical_epsilon_LOG_bicri,Z_LOG_bicri,Y_LOG_bicri] = corr_logdet_bicri(X_REL,C,k,epsilon,M,alpha);

%convex iteration for RCCP with alpha > 0
[X_CI_bicri,hodnost_CI_bicri,hodnost_final_CI_bicri,g_CI_bicri,g_final_CI_bicri,norma_CI_bicri,norm_final_CI_bicri,cas_CI_bicri,t_CI_bicri,s_CI_bicri,empirical_epsilon_CI_bicri,sum_eig_final_CI_bicri,sum_eig_real_CI_bicri,Z_CI_bicri,Y_CI_bicri] = corr_cvx_iter_bicri(C,k,epsilon,M,alpha);

%logdet for RCFP modified
[X_LOG_MOD,hodnost_LOG_MOD,hodnost_final_LOG_MOD,g_LOG_MOD,g_final_LOG_MOD,norma_LOG_MOD,norm_final_LOG_MOD,cas_LOG_MOD,t_LOG_MOD,s_LOG_MOD,empirical_epsilon_LOG_MOD,Z_LOG_MOD,Y_LOG_MOD] = corr_logdet_mod(X_REL,C,k,epsilon,M,gamma);

%convex iteration for RCFP modified
[X_CI_MOD,hodnost_CI_MOD,hodnost_final_CI_MOD,g_CI_MOD,g_final_CI_MOD,norma_CI_MOD,norm_final_CI_MOD,cas_CI_MOD,t_CI_MOD,s_CI_MOD,empirical_epsilon_CI_MOD,sum_eig_final_CI_MOD,sum_eig_real_CI_MOD,Z_CI_MOD,Y_CI_MOD] = corr_cvx_iter_mod(C,k,epsilon,M,gamma);

%optimal decomposition-based rank reduction algorithm
[X_RRO,hodnost_RRO,hodnost_final_RRO,g_RRO,g_final_RRO,norma_RRO,norm_final_RRO,cas_RRO,t_RRO,s_RRO,empirical_epsilon_RRO] = corr_reduction_opt(C,k,epsilon,M);

%feasible decomposition-based rank reduction algorithm
[X_RRF,hodnost_RRF,hodnost_final_RRF,g_RRF,g_final_RRF,norma_RRF,norm_final_RRF,cas_RRF,t_RRF,s_RRF,empirical_epsilon_RRF, koniec] = corr_reduction_feas(C,k,epsilon,M);

%bisective algorithm with logdet for RCFP
[X_hat,hodnost_hat,g_hat,norma_hat,gamma,hodnost,iter,cas,hodnost_final_LOG_MOD,g_final_LOG_MOD,norm_final_LOG_MOD,cas_LOG_MOD,t_LOG_MOD,s_LOG_MOD,empirical_epsilon_LOG_MOD] = bisection_logdet(g_final_REL,g_final_LOG,X_REL,X_LOG,C,k,epsilon,M);

%bisective algorithm with convex iteration for RCFP
[X_hat,hodnost_hat,g_hat,norma_hat,gamma,hodnost,iter,cas,hodnost_final_CI_MOD,g_final_CI_MOD,norm_final_CI_MOD,cas_CI_MOD,t_CI_MOD,s_CI_MOD,empirical_epsilon_CI_MOD] = bisection_CI(g_final_REL,g_final_CI,X_REL,X_CI,C,k,epsilon,M);

%bisective algorithm with logdet for RCCP with alpha > 0
[X_hat,hodnost_hat,g_hat,norma_hat,gamma,hodnost,iter,cas,hodnost_final_LOG_MOD,g_final_LOG_MOD,norm_final_LOG_MOD,cas_LOG_MOD,t_LOG_MOD,s_LOG_MOD,empirical_epsilon_LOG_MOD] = bisection_logdet(g_final_REL,g_final_LOG_bicri,X_REL,X_LOG_bicri,C,k,epsilon,M);

%bisective algorithm with convex iteration for RCCP with alpha > 0
[X_hat,hodnost_hat,g_hat,norma_hat,gamma,hodnost,iter,cas,hodnost_final_CI_MOD,g_final_CI_MOD,norm_final_CI_MOD,cas_CI_MOD,t_CI_MOD,s_CI_MOD,empirical_epsilon_CI_MOD] = bisection_CI(g_final_REL,g_final_CI_bicri,X_REL,X_CI_bicri,C,k,epsilon,M);
