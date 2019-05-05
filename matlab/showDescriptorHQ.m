function showDescriptorHQ(X,T,nid,n,V,refScale,h,w)
% Visualize a function on the mesh (with tessellation)
% The mesh is represented by
%   X - Vertex positions
%   T - Vertex indices on each triangle

if h*w+1 > size(V,2)
    fprintf('Error: number of eigenvectors not sufficient\n');
    return;
end

% Build refined mesh (the same resolution as n=6 by default)
nt = size(T,1);
res = 6;
[Xf,Tf,nidf] = meshDivide(X,T,res);
coeff = coordCoeff(res,n,refScale)*showBasis(n,refScale,0);
Vf = zeros(size(Xf,1),size(V,2));
for ti = 1:nt
    Vf(nidf(ti,:),:) = coeff*V(nid(ti,:),:);
end

% Draw refined mesh
figure;
cameratoolbar; cameratoolbar('SetCoordSys','none');
for ri = 1:h
    for ci = 1:w
        ind = (ri-1)*w+ci;
        subplot(h,w,ind);
        patch('vertices',Xf,'Faces',Tf,'FaceColor','interp','CData',Vf(:,ind+1),'edgecolor','none');
        % view(0,0);  % Needed when visualizing sphere meshes
        axis equal;
        axis off;
    end
end

end
