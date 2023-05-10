function[X,hodnost,hodnost_final,g,g_final,norma,norm_final,cas,t,s,empirical_epsilon] = corr_reduction_opt(C,k,epsilon,M)

n = size(C,1);

tic;
%initialize
[X,bla1,bla2,bla3,bla4,bla5,Z] = corr_relaxation(C,k,epsilon);
e = eye(3*n);

g = sum(diag(Z));
norma = norm(C-X,'fro');

%set counter of iteration
t = 0;

%set counter of the same iterations
s = 0;

Y = [X, zeros(n,2*n); zeros(n,n), X, C-X; zeros(n,n), (C-X)', Z];
hodnost = sum(eig(X)>epsilon);

while (hodnost(end) > k && s < M)
    
r = sum(eig(Y)>epsilon);

%rozklad matice X
[s,v] = svd(Y);
V = s*v^(1/2);
V = V(:,1:r);

pom = [zeros(n,n), ones(n,2*n); ones(2*n,n), zeros(2*n,2*n)];
C = [zeros(2*n,3*n); zeros(n,2*n), eye(n)];
%find direction matrix delta
cvx_begin
variable delta(r,r) symmetric
minimize sum(diag(delta))
for i = 1:n
    sum(diag(V'*e(:,i)*e(:,i)'*V*delta)) == 0;
end
pom.*(V*delta*V') == pom.*[eye(n), zeros(n,2*n); zeros(2*n,n), zeros(2*n,2*n)];
sum(diag(V'*C*V*delta)) == 0;
cvx_end

%maximum magnitude eigenvalue of delta
vh = eig(delta);
lambda = vh(find(abs(vh)==max(abs(vh))));

%step size
alpha = -1/lambda;

%counter of iterations
t = t + 1;

%shifted matrix
Y = V*(eye(r)+alpha*delta)*V';
X = Y(1:n,1:n);

%rank determination
hodnost = [hodnost;sum(eig(X)>epsilon)];
norma = [norma;norm(C-X,'fro')];

%pocitadlo iteracii, v ktorych nezmenena hodnosti
if hodnost(end) == hodnost(end-1)
    s = s + 1;
else
    s = 1;
end
end

%define outputs
cas = toc;
hodnost_final = hodnost(end);
g_final = sum(diag(Z));
norm_final = norm(C-X,'fro');
vh = eig(X);
empirical_epsilon = vh(n-k);
end