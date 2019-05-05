% Discrete Laplacian-Beltrami operator using higher-order FEM
% by Beichen Li
% Known limitation: calculation of basis functions is reliable only with
% n <= 10, due to numerical issues.

% Global parameters
n = 1;                      % Order of finite element
displayBasis = 0;           % Show a figure of all basis functions
tableInverse = 0;           % Store denominators in integral table
refScale = 1;               % Scale of reference triangle
meshPath = '../meshes';     % Path to mesh files
meshName = 'human_coarse_rot.off';  % Name of mesh file
mode = 'eval';              % Result processing mode ('eval' or 'vis')

% Evaluation ('eval') mode parameters
gtDataPath = '../gtdata';   % Path to ground-truth data folder
gtReadFile = 1;             % Read ground-truth data from external file
gtWriteFile = 0;            % Write ground-truth data to external file
eEig = 100;                 % Number of eigenfunctions
gtEig = 100;                % Number of ground-truth eigenfunctions
gtRes = 20;                 % Resolution of ground-truth refined mesh

% Visualization ('vis') mode parameters
visQuality = 1;             % Enable higher visualization quality when n <= 3
hEig = 4;                   % Number of eigenfunctions per column
wEig = 6;                   % Number of eigenfunctions per row

% Read input object
[X,T] = readOff(join([meshPath,'/',meshName]));

% Solve eigenvalues, eigenfunctions and obtain subdivided mesh
nEig = max(hEig*wEig+1,eEig);
[V,~,Xa,Ta,nid,~] = LaplaceBeltrami(X,T,n,nEig,refScale,displayBasis,tableInverse);

% Evaluate error in 'eval' mode
if strcmp(mode,'eval')
    % Build refined mesh
    [Xgt,Tgt,nidgt] = meshDivide(X,T,gtRes);
    nvgt = size(Xgt,1);
    
    % Calculate eigenfunction values on refined mesh
    nt = size(T,1);
    coeff = coordCoeff(gtRes,n,refScale)*showBasis(n,refScale,0);
    Ve = zeros(nvgt,eEig);
    for ti = 1:nt
        Ve(nidgt(ti,:),:) = coeff*V(nid(ti,:),1:eEig);
    end
    
    % Compute error
    gtDataFile = join([gtDataPath,'/',meshName(1:end-4),'_gt.mat']);
    err = evalError(Ve,gtEig,gtReadFile,gtWriteFile,gtDataFile);

% Show rendered mesh in 'vis' mode
elseif strcmp(mode,'vis')
    if visQuality > 0 && (n == 2 || n == 3)
        showDescriptorHQ(X,T,nid,n,V,refScale,hEig,wEig);
    else
        showDescriptor(Xa,Ta,V,hEig,wEig);
    end

% Undefined mode
else
    printf('Error: undefined mode');
end
