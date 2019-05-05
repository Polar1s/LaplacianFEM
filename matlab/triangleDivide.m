function Xt = triangleDivide(X,n)
% Generate the subdivided mesh of a triangle in certain resolution
%   X - 3*2 or 3*3 matrix of vertex coordinates

dim = size(X,2);
nv = (n+1)*(n+2)/2;

% Fill point positions
Xt = zeros(nv,dim);
Xt([1 n*(n+1)/2+1 nv],:) = X;
seq = zeros(3,n-1);
seq(1,:) = (1:n-1).*(2:n)./2+1;
seq(2,:) = (2:n)+(n*(n+1)/2);
seq(3,:) = seq(1,n-1:-1:1)+(n-1:-1:1);
for i = 1:3
    j = mod(i,3)+1;
    Xt(seq(i,:),:) = ((1:n-1)'*X(j,:)+(n-1:-1:1)'*X(i,:))/n;
end
p = X(2,:);
vb = X(1,:)-p;
va = X(3,:)-p;
for i = 1:n-2
    k = (i+1)*(i+2)/2;
    Xt(k+2:k+i+1,:) = p+((1:i)'*va+(n-1-i)*vb)/n;
end

end
