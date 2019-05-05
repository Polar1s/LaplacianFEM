function A = integralTable(n,refScale,denom)
% Table of int_Tri x^n*y^m dA (Tri = reference triangle)

N = 2*n+1;
A = ones(N);
P = ones(N);
for j = 1:N
    for i = 1:N+1-j
        if i > 1 && j > 1
            P(i,j) = P(i-1,j)+P(i,j-1);
        end
        A(i,j) = P(i,j)*(i+j-1)*(i+j)/refScale^(i+j);
    end
end
if denom == 0
    A = 1./A;
end
