function basis = showBasis(n,refScale,display)
% Calculate basis functions of higher-order elements

% Degree of freedom
dof = (n+1)*(n+2)/2;
basis = coordCoeff(n,n,refScale)\eye(dof); % Numerical issues with large n (n~20)

% Prepare triangular mesh
if display
    res = 20;
    Z = coordCoeff(res,n,refScale)*basis;
    X = triangleDivide([0.5 1;0 0;1 0]*refScale,res);
    T = zeros(n*n,3);
    T(1,:) = [1 2 3];
    for i = 2:n
        j = i*(i-2)+2;
        k = i*(i+1)/2;
        T(j:2:i*i,:) = [k-i+1:k;k+1:k+i;k+2:k+i+1]';
        T(j+1:2:i*i,:) = [k-i+1:k-1;k+2:k+i;k-i+2:k]';
    end

    % Draw mesh
    figure;
    ind = 1;
    for i = 1:n+1
        for j = 1:i
            trisurf(T,X(:,1)+(j-1)*refScale,X(:,2)+(n-i+1)*refScale,Z(:,ind));
            ind = ind+1;
            hold on;
        end
    end
    axis off;
    view(0,90);
    shading interp;
end

end
