function [E,Te] = edgeIndex(T)
% Index edges of input mesh and establish links between edges and triangles
nt = size(T,1);

% Count number of edges
ne = 0;
for j = 1:3
    k = mod(j,3)+1;
    for i = 1:nt
        if T(i,j) < T(i,k)
            ne = ne+1;
        end
    end
end

% Collect and sort edges and perform duplication check
E = zeros(ne,2);
ei = 1;
for j = 1:3
    k = mod(j,3)+1;
    for i = 1:nt
        if T(i,j) < T(i,k)
            E(ei,:) = [T(i,j) T(i,k)];
            ei = ei+1;
        end
    end
end
E(:) = sortrows(E);
for i = 1:ne-1
    if isequal(E(i,:),E(i+1,:))
        printf('Error: duplicate edges');
    end
end

% Link triangles to edges
Te = zeros(nt,3);
for j = 1:3
    k = mod(j,3)+1;
    for i = 1:nt
        u = min(T(i,j),T(i,k));
        v = max(T(i,j),T(i,k));
        ei = searchEdge(E,u,v);
        if ei == 0
            printf('Error: edge (%d,%d) not found in list',u,v);
        else
            Te(i,j) = ei;
        end
    end
end
end

function ei = searchEdge(E,u,v)
% Search for an edge (u,v) (u <= v) in list E
ne = size(E,1);
l = 1;
r = ne;
ei = 0;
while l <= r
    m = floor((l+r)/2);  % Do not use idivide since it is very slow
    if E(m,1) == u && E(m,2) == v
        ei = m;
        break;
    elseif E(m,1) > u || E(m,1) == u && E(m,2) > v
        r = m-1;
    else
        l = m+1;
    end
end
end
