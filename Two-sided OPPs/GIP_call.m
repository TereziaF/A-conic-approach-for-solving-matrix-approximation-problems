%Graph Isomorphism Problem
C = [0 1 0 0 1 0;
     1 0 1 0 0 0;
     0 1 0 1 0 1;
     0 0 1 0 1 0;
     1 0 0 1 0 1;
     0 0 1 0 1 0];
 
A = [0 0 0 1 0 1;
     0 0 0 1 0 1;
     0 0 0 1 1 0;
     1 1 1 0 0 0;
     0 0 1 0 0 1;
     1 1 0 0 1 0];
    
%define inputs for methods
epsilon = 10^(-6);
alpha = 100;
M = 20;
m = size(A,2);
n = size(C,1);
p = n;
q = m;

%relaxation
tic;
[X0,t0,Y0,V0] = SDPrelax_L1_PERMUTATION_GIP(A, C, m, n, p, q) 
cas = toc;


ucel = norm(X0*C-A*X0,1);                    %optimal value
orth = norm(X0'*X0-eye(m),'fro');            %orthogonality criterion
sum_row = norm(X0*ones(n)-ones(m),1);        %test if sum of each row equals 1
sum_col = norm(X0'*ones(m)-ones(n),1);       %test if sum of each column equals 1
max5 = sort(vec(X0));                        %ordering elements
o = max5(end-(n-1):end);                     %n largest elements
max_dev(opakuj) = norm(o-ones(n),1);         %test if n largest elements equal 1
z = max5(1:end-n);                           %n(n-1) smallest elements
min_dev(opakuj) = norm(z-zeros(n^2-n),1);    %test if n(n-1) smallest elements equal 0
X0 = full(X0); 
hodnost = sum(eig(V0)>epsilon);              %rank
vh = eig(V0);                                %eigenvalues
emp_eps = vh(end-n+1);                       %empirical epsilon                 

%modified logdet heuristic
if (hodnost > n)
gamma = epsilon;
tic;
[X,hodnost,hodnost_final,g,g_final,norma,norm_final,cas,t,s,empirical_epsilon,t_var,S,V] = OPP_logdet_mod_PERMUTATION_GIP(X0,V0,C,A,k,epsilon,M,gamma)
Lcas = toc;
end

Lucel = norm(X*C-A*X,1);                      %optimal value
Lorth = norm(X'*X-eye(m),'fro');              %orthogonality criterion
Lsum_row = norm(X*ones(n)-ones(m),1);         %test if sum of each row equals 1
Lsum_col = norm(X'*ones(m)-ones(n),1);        %test if sum of each column equals 1
max5 = sort(vec(X));                          %ordering elements
o = max5(end-(n-1):end);                      %n largest elements
Lmax_dev = norm(o-ones(n),1);                 %test if n largest elements equal 1
z = max5(1:end-n);                            %n(n-1) smallest elements
Lmin_dev = norm(z-zeros(n^2-n),1);            %test if n(n-1) smallest elements equal 0
X = full(X); 
Lhodnost = sum(eig(V)>epsilon);               %rank

%modified convex iteration
tic;
[X_CI,hodnost_CI,hodnost_final_CI,g_CI,g_final_CI,norma_CI,norm_final_CI,cas_CI,t_CI,s_CI,empirical_epsilon_CI,sum_eig_final,sum_eig_real,t_var_CI,S_CI,V_CI] = OPP_cvx_iter_mod_L1_PERMUTATION_GIP(C,A,k,epsilon,M,gamma);
Ccas = toc;

Cucel = norm(X_CI*C-A*X_CI,1);                %optimal value
Corth = norm(X_CI'*X_CI-eye(m),'fro');        %orthogonality criterion
Csum_row = norm(X_CI*ones(n)-ones(m),1);      %test if sum of each row equals 1
Csum_col = norm(X_CI'*ones(m)-ones(n),1);     %test if sum of each column equals 1
max5 = sort(vec(X_CI));                       %ordering elements
o = max5(end-(n-1):end);                      %n largest elements
Cmax_dev = norm(o-ones(n),1);                 %test if n largest elements equal 1
z = max5(1:end-n);                            %n(n-1) smallest elements
Cmin_dev = norm(z-zeros(n^2-n),1);            %test if n(n-1) smallest elements equal 0
X_CI = full(X_CI);                            
Chodnost = sum(eig(V_CI)>epsilon);            %rank