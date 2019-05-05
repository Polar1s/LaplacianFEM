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
err = zeros(n,eEig);
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

% Evaluate speed at each order from 1 to n
tc = zeros(n,1);
for k = 1:n
    fprintf('Evaluating speed at n = %d ...\n',k);
    f = @() LaplaceBeltrami(X,T,k,eEig,refScale,0,tableInverse);
    tc(k) = timeit(f);
end

end
