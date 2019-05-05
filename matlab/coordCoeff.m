function C = coordCoeff(res,n,refScale)
% Calculate coefficient matrix for solving basis functions and computing
% basis function values at points

nvt = (res+1)*(res+2)/2;
dof = (n+1)*(n+2)/2;
C = ones(nvt,dof);
x = zeros(nvt,1);
y = zeros(nvt,1);
for i = 1:res+1
    x((i-1)*i/2+1:i*(i+1)/2) = (0:i-1)*refScale/res;
    y((i-1)*i/2+1:i*(i+1)/2) = (1-(i-1)/res)*refScale;
end
for i = 2:n+1
    k = i*(i+1)/2;
    C(:,k-i+1:k) = (x.^(i-1:-1:0)).*(y.^(0:i-1));
end

end
