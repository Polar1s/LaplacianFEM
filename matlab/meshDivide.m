function [Xa,Ta,nid] = meshDivide(X,T,n)
% Mesh subdivision for visualization

nv = size(X,1);
nt = size(T,1);
dof = (n+1)*(n+2)/2;
[E,Te] = edgeIndex(T);

% Build node index matrix
nid = zeros(nt,dof);
ne = size(E,1);
nve = nv+ne*(n-1);
nTot = nve+nt*(n-2)*(n-1)/2;
seq = zeros(3,n-1);
seq(1,:) = (1:n-1).*(2:n)./2+1;
seq(2,:) = (2:n)+(n*(n+1)/2);
seq(3,:) = seq(1,n-1:-1:1)+(n-1:-1:1);

for ti = 1:nt
    nid(ti,[1 n*(n+1)/2+1 dof]) = T(ti,:);
    for i = 1:3
        j = mod(i,3)+1;
        if T(ti,i) < T(ti,j)
            nid(ti,seq(i,:)) = (1:n-1)+(nv+(Te(ti,i)-1)*(n-1));
        else
            nid(ti,seq(i,:)) = (n-1:-1:1)+(nv+(Te(ti,i)-1)*(n-1));
        end
    end
    for i = 1:n-2
        k = (i+1)*(i+2)/2;
        nid(ti,k+2:k+i+1) = (k-2*i:k-i-1)+(nve+(ti-1)*(n-2)*(n-1)/2);
    end
end

% Augment vertex and triangle set
Xa = zeros(nTot,3);
Ta = zeros(nt*n*n,3);
Tt = zeros(n*n,3);
Tt(1,:) = [1 2 3];
for i = 2:n
    j = i*(i-2)+2;
    k = i*(i+1)/2;
    Tt(j:2:i*i,:) = [k-i+1:k;k+1:k+i;k+2:k+i+1]';
    Tt(j+1:2:i*i,:) = [k-i+1:k-1;k+2:k+i;k-i+2:k]';
end
for ti = 1:nt
    Xt = triangleDivide(X(T(ti,:),:),n);
    Xa(nid(ti,:),:) = Xt;
    off = (ti-1)*n*n;
    Ta(off+1:off+n*n,:) = reshape(nid(ti,Tt),size(Tt));
end

end
