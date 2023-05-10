%approximation correlation matrix generation - diagonal 1, off-diagonal in
%[-1,1]

function[C] = approx_corr(n)
C = ones(n)-2*rand(n);
for i = 1:n
    for j = 1:n
        if i < j
            C(i,j) = C(j,i);
        elseif i == j
                C(i,i) = 1;
        end
    end
end
end