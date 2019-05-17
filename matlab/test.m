% Discrete Laplacian-Beltrami operator using higher-order FEM
% by Beichen Li
% Known limitation: calculation of basis functions is reliable only with
% n <= 10, due to numerical issues.

% Global parameters
n = 5;                     % Order of finite element
displayBasis = 0;           % Show a figure of all basis functions
tableInverse = 0;           % Store denominators in integral table
refScale = 1;               % Scale of reference triangle
meshPath = '../meshes';     % Path to mesh files
meshName = '166_lite.off';  % Name of mesh file
mode = 'vis2';              % Result processing mode ('eval' or 'vis')

% Evaluation ('eval') mode parameters
gtDataPath = '../gtdata';   % Path to ground-truth data folder
gtReadFile = 0;             % Read ground-truth data from external file
gtWriteFile = 1;            % Write ground-truth data to external file
eEig = 100;                 % Number of eigenfunctions
gtEig = 300;                % Number of ground-truth eigenfunctions
gtRes = 20;                 % Resolution of ground-truth refined mesh

% Visualization ('vis') mode parameters
visQuality = 1;             % Enable higher visualization quality when n <= 3
hEig = 4;                   % Number of eigenfunctions per column
wEig = 6;                   % Number of eigenfunctions per row

% Read input object
[X,T] = readOff(join([meshPath,'/',meshName]));

% Evaluate error in 'eval' mode
if strcmp(mode,'eval')
    gtDataFile = join([gtDataPath,'/',meshName(1:end-4),'_gt.mat']);
    [tc,err] = evalError(X,T,n,eEig,gtRes,gtEig,refScale,gtReadFile,gtWriteFile,gtDataFile,tableInverse);
    
    resFile = join(['../results/',meshName(1:end-4),'_res.mat']);
    save(resFile,'tc','err');

% Show eigenfunctions in 'vis' mode
elseif strcmp(mode,'vis1')
    nEig = hEig*wEig+1;
    [V,~,Xa,Ta,nid,~] = LaplaceBeltrami(X,T,n,nEig,refScale,displayBasis,tableInverse);
    
    if visQuality > 0 && (n == 2 || n == 3)
        showDescriptorHQ(X,T,nid,n,V,refScale,hEig,wEig);
    else
        showDescriptor(Xa,Ta,V,hEig,wEig);
    end
elseif strcmp(mode,'vis2')
    [V,~,Xa,Ta,nid,~] = LaplaceBeltrami(X,T,n,25,refScale,displayBasis,tableInverse);
    inds = [1 5 8 13 25];
    
    if visQuality > 0 && (n == 2 || n == 3)
        showDescriptorHQ(X,T,nid,n,V(:,inds),refScale,1,4);
    else
        showDescriptor(Xa,Ta,V(:,inds),1,4);
    end

% Undefined mode
else
    printf('Error: undefined mode');
end
