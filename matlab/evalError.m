function err = evalError(Ve,gtEig,readFile,writeFile,fileName)
% Evaluate the error of a set of eigenfunctions on ground-truth mesh

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

% Compute error
Verr = Vgt*(Vgt'*Mgt*Ve)-Ve;
err = sqrt(diag(Verr'*Mgt*Verr));

end
