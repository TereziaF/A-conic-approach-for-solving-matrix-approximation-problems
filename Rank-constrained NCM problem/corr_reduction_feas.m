function[X,hodnost,hodnost_final,g,g_final,norma,norm_final,cas,t,s,empirical_epsilon,koniec,g_final_SDP,norm_final_SDP, hodnost_SDP,cas_SDP] = corr_reduction_feas(C,k,epsilon,M)

n = size(C,1);
koniec = 0;

tic;
%initialize
[X,hodnost_SDP,g_final_SDP,norm_final_SDP,cas_SDP,bla5,Z] = corr_relaxation(C,k,epsilon);
e = eye(n);
vh = eig(X);
hodnost = sum(eig(X)>epsilon);
g = sum(diag(Z));
norma = norm(C-X,'fro');
emp = vh(n-k);

%set counter of iteration
t = 0;

%set counter of the same iterations
s = 0;

while (hodnost(end) > k && s < M)
    
r = hodnost(end);

%rozklad matice X
[s,v] = svd(X);
V = s*v^(1/2);
V = V(:,1:r);

%find direction matrix delta
cvx_begin
variable delta(r,r) symmetric
minimize sum(diag(delta))
for i = 1:n
    sum(diag(V'*e(:,i)*e(:,i)'*V*delta)) == 0;
end
sum(diag(delta)) >= 0;
cvx_end

if isnan(cvx_optval)
    cas = toc;
    hodnost_final = hodnost(end);
    g_final = g(end);
    norm_final = norma(end);
    empirical_epsilon = emp(end);
    koniec = 1;
    break;
else
%maximum magnitude eigenvalue of delta
vh = eig(delta);
lambda = vh(find(abs(vh)==max(abs(vh))));

%step size
alpha = -1/lambda;

%counter of iterations
t = t + 1;

%shifted matrix
X = V*(eye(r)+alpha*delta)*V';

%rank determination
hodnost = [hodnost;sum(eig(X)>epsilon)];
norma = [norma;norm(C-X,'fro')];
g = [g;sum(diag(Z))];
vh = eig(X);
emp = vh(n-k);

%pocitadlo iteracii, v ktorych nezmenena hodnosti
if hodnost(end) == hodnost(end-1)
    s = s + 1;
else
    s = 1;
end


end
end

%define outputs
cas = toc;
hodnost_final = hodnost(end);
g_final = g(end);
norm_final = norma(end);
empirical_epsilon = emp(end);
end