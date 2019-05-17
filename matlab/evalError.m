function [tc,err] = evalError(X,T,n,eEig,gtRes,gtEig,refScale,readFile,writeFile,fileName,tableInverse)
% Evaluate the error of a set of eigenfunctions on ground-truth mesh

% Build refined mesh
[Xgt,Tgt,nidgt] = meshDivide(X,T,gtRes);
nvgt = size(Xgt,1);

% Compute ground-truth eigenfunction space or read from file
if readFile > 0
    load(fileName,'Vgt','Mgt');
else
    [Vgt,~,~,~,~,Mgt] = LaplaceBeltrami(Xgt,Tgt,1,gtEig,1,0,0);
end
% Save ground-truth data
if writeFile > 0
    save(fileName,'Vgt','Mgt');
end

% Evaluate error at each order from 1 to n
err = zeros(2*n,eEig);
for k = 1:n
    fprintf('Evaluating error at n = %d ...\n',k);
    
    % Calculate eigenfunction values on refined mesh
    [V,~,~,~,nid,~] = LaplaceBeltrami(X,T,k,eEig,refScale,0,tableInverse);
    
    coeff = coordCoeff(gtRes,k,refScale)*showBasis(k,refScale,0);
    Ve = zeros(nvgt,eEig);
    for ti = 1:size(T,1)
        Ve(nidgt(ti,:),:) = coeff*V(nid(ti,:),1:eEig);
    end

    % Compute error
    Verr = Vgt*(Vgt'*Mgt*Ve)-Ve;
    err(k,:) = sqrt(diag(Verr'*Mgt*Verr));
end
for k = 1:n
    fprintf('Evaluating error at res = %d ...\n',k);
    
    % Calculate refined mesh
    [Xr,Tr,nidr] = meshDivide(X,T,k);
    [Vr,~,~,~,~,~] = LaplaceBeltrami(Xr,Tr,1,eEig,refScale,0,tableInverse);
    
    P = interpMatrix(gtRes,k);
    Ve = zeros(nvgt,eEig);
    for ti = 1:size(T,1)
        Ve(nidgt(ti,:),:) = P*Vr(nidr(ti,:),1:eEig);
    end
    
    % Compute error
    Verr = Vgt*(Vgt'*Mgt*Ve)-Ve;
    err(n+k,:) = sqrt(diag(Verr'*Mgt*Verr));
end

% Evaluate speed at each order from 1 to n
tc = zeros(2*n,1);
for k = 1:n
    fprintf('Evaluating speed at n = %d ...\n',k);
    f = @() LaplaceBeltrami(X,T,k,eEig,refScale,0,tableInverse);
    tc(k) = timeit(f);
end
for k = 1:n
    fprintf('Evaluating speed at res = %d ...\n',k);
    f = @() LBRefineTest(X,T,k,eEig,refScale,tableInverse);
    tc(n+k) = timeit(f);
end

end

function P = interpMatrix(nt,ns)
% Interpolation matrix for calculating values of refined linear FEM
dofs = (ns+1)*(ns+2)/2;
doft = (nt+1)*(nt+2)/2;
P = zeros(doft,dofs);

% Fill in entries
for i = 0:nt
    v = ns-i*ns/nt;
    fv = min(ns-1,floor(v));
    rv = v-fv;
    for j = 0:i
        u = j*ns/nt;
        fu = min(ns-fv-1, floor(u));
        ru = u-fu;
        id = (ns-fv-1)*(ns-fv)/2+fu+1;
        idt = i*(i+1)/2+j+1;
        if ru+rv <= 1  % Lower triangle
            P(idt,[id id+ns-fv id+ns-fv+1]) = [rv 1-ru-rv ru];
        elseif ru <= 1 && rv <= 1  % Upper triangle
            P(idt,[id id+1 id+ns-fv+1]) = [1-ru ru+rv-1 1-rv];
        else  % Error
            fprintf('Error: interpolation out of boundary!\n');
        end
    end
end

end

function [V,d,Xa,Ta,nid,M] = LBRefineTest(X,T,res,nEig,refScale,tableInverse)
% Wrapped function for speed test
[Xa,Ta,nid] = meshDivide(X,T,res);
[V,d,~,~,~,M] = LaplaceBeltrami(Xa,Ta,1,nEig,refScale,0,tableInverse);
end
