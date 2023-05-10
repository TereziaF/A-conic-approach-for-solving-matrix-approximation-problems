tabulka = [...
    "orth_final_REL","hodnost_final_REL","g_final_REL","norm_final_REL","cas_REL","emp_eps_REL",...
   "orth_final_LOG_bicri","hodnost_final_LOG_bicri","g_final_LOG_bicri","norm_final_LOG_bicri","cas_LOG_bicri","emp_eps_LOG_bicri",...
   "orth_final_CI_bicri","hodnost_final_CI_bicri","g_final_CI_bicri","norm_final_CI_bicri","cas_CI_bicri","empirical_epsilon_CI_bicri",...
    "orth_final_LOG_MOD","hodnost_final_LOG_MOD","g_final_LOG_MOD","norm_final_LOG_MOD","cas_LOG_MOD","empirical_epsilon_LOG_MOD",...
   "orth_final_CI_MOD","hodnost_final_CI_MOD","g_final_CI_MOD","norm_final_CI_MOD","cas_CI_MOD","empirical_epsilon_CI_MOD",...
    "orth_hat_LOG_bisection", "hodnost_hat_LOG","g_hat_LOG","norma_hat_LOG", "iter_LOG", "cas_BI_LOG",...
   "orth_hat_CI_bisection","hodnost_hat_CI","g_hat_CI","norma_hat_CI", "iter_CI", "cas_BI_CI"];
m = 4;
n = 4;
p = 10;
q = 3;
epsilon = 10^(-6);
tol = epsilon;
M = 10;
maxiter = M;

 k=m;
 
%generating 100 problems
for opakuj = 1:100
rng(opakuj)
A = randn(p,m);
B = randn(n,q);

%generating oblique solution
X_sol = randn(m,n);
X_sol = X_sol*diag(diag(X_sol'*X_sol).^(-1/2))

%delta = 0;                    %zero gap
delta = randn(p,q);            %non-zero gap
C = A*X_sol*B + 0.5*delta;

W = ones(p,q);

tic;
[X_rel,t_rel,S_rel,V_rel,G_rel] = SDPrelax_OBLIQUE_Linfty(A, B, C, W, m, n, p, q);
cas_REL = toc;                                              %time
hodnost_final_REL = sum(eig(V_rel)>epsilon);                %rank
g_final_REL = t_rel;                                        %optimal value of reformulation
norm_final_REL = norm(C-A*X_rel*B,Inf);                     %optimal value
orth_final_REL = norm(diag(X_rel'*X_rel)-ones(n,1),1);      %oblique criterion
vh = eig(V_rel);                                            %eigenvalues
emp_eps_REL = vh(n+1);                                      %empirical epsilon

if (hodnost_final_REL > m)
%logdet bicri
alpha = 10;
[X_LOG_bicri,hodnost_LOG_bicri,hodnost_final_LOG_bicri,g_LOG_bicri,g_final_LOG_bicri,norma_LOG_bicri,norm_final_LOG_bicri,cas_LOG_bicri,t_LOG_bicri,s_LOG_bicri,empirical_epsilon_LOG_bicri,t_var_LOG_bicri,S_LOG_bicri,V_LOG_bicri,G_LOG_bicri] = logdet_bicri_Linfty_OBLIQUE(X_rel,V_rel,C,A,B,W,m,epsilon,M,alpha);
orth_final_LOG_bicri = norm(diag(X_LOG_bicri'*X_LOG_bicri)-ones(n,1),1);   %oblique criterion
hodnost_final_LOG_bicri                                                    %rank

%convex iter. bicri
alpha = 10;
[X_CI_bicri,hodnost_CI_bicri,hodnost_final_CI_bicri,g_CI_bicri,g_final_CI_bicri,norma_CI_bicri,norm_final_CI_bicri,cas_CI_bicri,t_CI_bicri,s_CI_bicri,empirical_epsilon_CI_bicri,sum_eig_final_CI_bicri,sum_eig_real_CI_bicri,t_var_CI_bicri,S_CI_bicri,V_CI_bicri,G_CI_bicri] = cvx_iter_bicri_Linfty_OBLIQUE(C,A,B,W,m,epsilon,M,alpha);
orth_final_CI_bicri = norm(diag(X_CI_bicri'*X_CI_bicri)-ones(n,1),1);      %oblique criterion
hodnost_final_CI_bicri                                                     %rank

if (hodnost_final_LOG_bicri == m && abs(t_rel-g_final_LOG_bicri)>epsilon)
    %bisective algorithm with logdet for RCCP with alpha > 0
[X_hat,hodnost_hat_LOG,g_hat_LOG,norma_hat_LOG,gamma,hodnost,iter_LOG,cas,hodnost_final_LOG_MOD,g_final_LOG_MOD,norm_final_LOG_MOD,cas_LOG_MOD,t_LOG_MOD,s_LOG_MOD,empirical_epsilon_LOG_MOD] = bisection_logdet_Linfty_OBLIQUE(t_rel,g_final_LOG_bicri,X_rel,X_LOG_bicri,V_rel,V_LOG_bicri,C,A,B,W,m,epsilon,M);
orth_hat_LOG_bisection = norm(diag(X_hat'*X_hat)-ones(n,1),1);             %oblique criterion
cas_BI_LOG = mean(cas_LOG_MOD(2:end))*iter_LOG;                            %time
hodnost_hat_LOG                                                            %rank
end

if (hodnost_final_CI_bicri == m && abs(t_rel-g_final_CI_bicri)>epsilon)
    %bisective algorithm with convex iteration for RCCP with alpha > 0
[X_hat,hodnost_hat_CI,g_hat_CI,norma_hat_CI,gamma,hodnost,iter_CI,cas,hodnost_final_CI_MOD,g_final_CI_MOD,norm_final_CI_MOD,cas_CI_MOD,t_CI_MOD,s_CI_MOD,empirical_epsilon_CI_MOD] = bisection_CI_Linfty_OBLIQUE(t_rel,g_final_CI_bicri,X_rel,X_CI_bicri,V_rel,V_CI_bicri,C,A,B,W,m,epsilon,M);
orth_hat_CI_bisection = norm(diag(X_hat'*X_hat)-ones(n,1),1);              %oblique criterion
cas_BI_CI = mean(cas_CI_MOD(2:end))*iter_CI;                               %time
hodnost_hat_CI                                                             %rank
end


%logdet for RCFP modified
gamma = t_rel;
[X_LOG_MOD,hodnost_LOG_MOD,hodnost_final_LOG_MOD,g_LOG_MOD,g_final_LOG_MOD,norma_LOG_MOD,norm_final_LOG_MOD,cas_LOG_MOD,t_LOG_MOD,s_LOG_MOD,empirical_epsilon_LOG_MOD,t_var,S,V,G] = logdet_mod_Linfty_OBLIQUE(X_rel,V_rel,C,A,B,W,m,tol,maxiter,gamma);
orth_final_LOG_MOD = norm(diag(X_LOG_MOD'*X_LOG_MOD)-ones(n,1),1);         %oblique criterion
hodnost_final_LOG_MOD                                                      %rank

%convex iteration for RCFP modified
gamma = t_rel;
[X_CI_MOD,hodnost_CI_MOD,hodnost_final_CI_MOD,g_CI_MOD,g_final_CI_MOD,norma_CI_MOD,norm_final_CI_MOD,cas_CI_MOD,t_CI_MOD,s_CI_MOD,empirical_epsilon_CI_MOD,sum_eig_final_CI_MOD,sum_eig_real_CI_MOD,t_var,S,V,G] = cvx_iter_mod_Linfty_OBLIQUE(C,A,B,W,m,tol,maxiter,gamma);
orth_final_CI_MOD = norm(diag(X_CI_MOD'*X_CI_MOD)-ones(n,1),1);            %oblique criterion  
hodnost_final_CI_MOD                                                       %rank

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tabulka = [tabulka;...
    orth_final_REL,hodnost_final_REL,g_final_REL,norm_final_REL,cas_REL,emp_eps_REL,...
   orth_final_LOG_bicri,hodnost_final_LOG_bicri,g_final_LOG_bicri,norm_final_LOG_bicri,cas_LOG_bicri,empirical_epsilon_LOG_bicri,...
   orth_final_CI_bicri,hodnost_final_CI_bicri,g_final_CI_bicri,norm_final_CI_bicri,cas_CI_bicri,empirical_epsilon_CI_bicri,...
   orth_final_LOG_MOD,hodnost_final_LOG_MOD,g_final_LOG_MOD,norm_final_LOG_MOD,sum(cas_LOG_MOD),mean(empirical_epsilon_LOG_MOD),...
   orth_final_CI_MOD,hodnost_final_CI_MOD,g_final_CI_MOD,norm_final_CI_MOD,sum(cas_CI_MOD),mean(empirical_epsilon_CI_MOD),...
    orth_hat_LOG_bisection, hodnost_hat_LOG,g_hat_LOG,norma_hat_LOG, iter_LOG, cas_BI_LOG,...
   orth_hat_CI_bisection,hodnost_hat_CI,g_hat_CI,norma_hat_CI, iter_CI, cas_BI_CI];
filename = ['OBLIQUE_Linfty_weighted_not0_P10N4M4Q3riadok' num2str(opakuj) '.csv'];
writematrix(tabulka,filename)
end