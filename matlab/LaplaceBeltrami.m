function [V,d] = LaplaceBeltrami(X,T,nid,basis,refScale,n,nTot,nEig,tableInverse)
% Solver for eigenvalues and eigenfunctions of Discrete Laplacian-Beltrami
% operator on a given mesh

nt = size(T,1);
dof = size(basis,1);
if dof ~= (n+1)*(n+2)/2
    printf('Error: input mesh has incorrect order');
    return;
end

% Constant matrices
IT = integralTable(n,refScale,tableInverse);
H = intProdMatrix(n,IT);
dCoeff = derivativeCoeff(n);

% Mass matrix per element (inner product of basis functions)
Me = intProd(basis,basis,n,H,tableInverse);

% Laplacian matrix base per element (inner product of gradient functions)
Le = zeros(4,dof,dof);
dBasis = derivative(basis,n,dCoeff);
Le(1,:,:) = intProd(dBasis(:,:,1),dBasis(:,:,1),n,H,tableInverse);
Le(2,:,:) = intProd(dBasis(:,:,1),dBasis(:,:,2),n,H,tableInverse);
Le(3,:,:) = intProd(dBasis(:,:,2),dBasis(:,:,1),n,H,tableInverse);
Le(4,:,:) = intProd(dBasis(:,:,2),dBasis(:,:,2),n,H,tableInverse);

% Construct mass matrix and Laplacian matrix
totData = nt*dof*dof;
Li = zeros(totData,1);
Lj = zeros(totData,1);
Lv = zeros(totData,1);
Mv = zeros(totData,1);
repMat = repmat(1:dof,[dof,1]);
repMatT = repMat';

for ti = 1:nt
    % Compute affine transform
    va = X(T(ti,3),:)-X(T(ti,2),:);
    vb = X(T(ti,1),:)-X(T(ti,2),:);
    a = norm(va,2);
    b = norm(vb,2);
    costh = dot(va,vb)/(a*b);
    sinth = sqrt(1-costh*costh);
    A = [1/a -costh/(a*sinth);0 1/(b*sinth)];
    
    % Fill matrix entries
    Li((1:dof*dof)+(ti-1)*dof*dof) = nid(ti,repMatT);
    Lj((1:dof*dof)+(ti-1)*dof*dof) = nid(ti,repMat);
    Lv((1:dof*dof)+(ti-1)*dof*dof) = reshape(sum(Le.*reshape(A*A',[4,1]),1),[dof,dof]);
    Mv((1:dof*dof)+(ti-1)*dof*dof) = Me*(a*b*sinth/refScale^2);
end

L = sparse(Li,Lj,Lv,nTot,nTot);
M = sparse(Li,Lj,Mv,nTot,nTot);

% Solve eigenfunctions of discrete Laplacian-Beltrami operator
[V,D] = eigs(L,M,nEig,'smallestabs');
d = diag(D);

% Normalize eigenfunctions
V(:) = real(V);
for i = 1:nEig
    V(:,i) = V(:,i)*sign(V(1,i));
end

end

%% Helper functions

function res = intProd(P,Q,n,H,tableInverse)
% Calculate the integral of product P(x)*Q(x) over reference triangle
% Both P(x) and Q(x) are nth-order polynomials represented in coefficients
dof = size(P,1);
if dof ~= (n+1)*(n+2)/2 || dof ~= size(Q,1)
    printf('Error: Input polynomial has incorrect order');
    res = 0;
else
    m = size(Q,2);
    HQ = zeros(dof,m);
    if tableInverse
        for i = 1:m
            HQ(:,i) = sum(Q(:,i)./H,2);
        end
    else
        HQ(:,:) = H*Q;
    end
    res = P'*HQ;
end
end

function H = intProdMatrix(n,T)
% Compute coefficient matrix needed for intProd function
dof = (n+1)*(n+2)/2;
H = zeros(dof);
ni = zeros(dof,1);
mi = zeros(dof,1);
for i = 1:n+1
    k = i*(i+1)/2;
    ni(k-i+1:k) = i:-1:1;
    mi(k-i+1:k) = 1:i;
end
for j = 1:dof
    for i = 1:dof
        H(i,j) = T(ni(i)+ni(j)-1,mi(i)+mi(j)-1);
    end
end
end

function Pd = derivative(P,n,C)
% Compute partial derivatives of nth-order polynomial function P(x)
dof = size(P,1);
m = size(P,2);
Pd = zeros(dof,m,2);
if dof ~= (n+1)*(n+2)/2
    printf('Error: Input polynomial has incorrect order');
else
    for i = 1:n
        k = i*(i+1)/2;
        Pd(k-i+1:k,:,1) = P(k+1:k+i,:);
        Pd(k-i+1:k,:,2) = P(k+2:k+i+1,:);
    end
    Pd(:,:,1) = Pd(:,:,1).*C(:,1);
    Pd(:,:,2) = Pd(:,:,2).*C(:,2);
end
end

function C = derivativeCoeff(n)
% Compute coefficient matrix used for partial derivatives
dof = (n+1)*(n+2)/2;
C = zeros(dof,2);
for i = 1:n
    k = i*(i+1)/2;
    C(k-i+1:k,:) = [i:-1:1;1:i]';
end
end