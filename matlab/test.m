% Discrete Laplacian-Beltrami operator using higher-order FEM
% by Beichen Li
% Known limitation: calculation of basis functions is reliable only with
% n <= 10, due to numerical issues.

% Global parameters
n = 2;                      % Order of finite element
displayBasis = 0;           % Show a figure of all basis functions
tableInverse = 0;           % Store denominators in integral table
refScale = 1;               % Scale of reference triangle
meshPath = '../meshes';     % Path to mesh files
meshName = 'human_coarse_rot.off';  % Name of mesh file
mode = 'vis';               % Result processing mode ('eval' or 'vis')

% Evaluation ('eval') mode parameters
eEig = 100;                 % Number of eigenfunctions
gtEig = 100;                % Number of ground-truth eigenfunctions

% Visualization ('vis') mode parameters
visQuality = 1;             % Enable higher visualization quality when n <= 3
hEig = 4;                   % Number of eigenfunctions per column
wEig = 6;                   % Number of eigenfunctions per row

% Read input object
[X,T] = readOff(join([meshPath,'/',meshName]));
[Xa,Ta,nid] = meshDivide(X,T,n);

% Main solver
basis = showBasis(n,refScale,displayBasis);
nEig = hEig*wEig+1;
nTot = size(Xa,1);
[V,~] = LaplaceBeltrami(X,T,nid,basis,refScale,n,nTot,nEig,tableInverse);

% Show rendered mesh
if visQuality > 0 && (n == 2 || n == 3)
    showDescriptorHQ(X,T,nid,basis,refScale,n,V,hEig,wEig);
else
    showDescriptor(Xa,Ta,V,hEig,wEig);
end
