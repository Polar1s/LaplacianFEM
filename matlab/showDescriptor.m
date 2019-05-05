function showDescriptor(X,T,V,h,w)
% Visualize a function on the mesh
% The mesh is represented by
%   X - Vertex positions
%   T - Vertex indices on each triangle

if h*w+1 > size(V,2)
    fprintf('Error: number of eigenvectors not sufficient\n');
    return;
end

% Draw augmented mesh
figure;
cameratoolbar; cameratoolbar('SetCoordSys','none');
for i = 1:h
    for j = 1:w
        subplot(h,w,(i-1)*w+j);
        patch('vertices',X,'Faces',T,'FaceColor','interp','CData',V(:,(i-1)*w+j+1),'edgecolor','none'); 
        % trimesh(T,X(:,1),X(:,2),X(:,3),V(:,(i-1)*w+j+1),'FaceColor','white','EdgeColor','interp');
        % view(0,0);  % Needed when visualizing sphere meshes
        axis equal;
        axis off;
    end
end
end
