    
currentFolder = pwd;

pathInput          = fullfile(currentFolder,'InputData');
pathAssembly       = fullfile(currentFolder,'Assembly');
pathErrors         = fullfile(currentFolder,'Errors');
pathMeshGeneration = fullfile(currentFolder,'MeshErrorAnalysis');
pathFEspace        = fullfile(currentFolder,'FEspace');
pathPostProcessing = fullfile(currentFolder,'PostProcessing');
pathMultiGrid      = fullfile(currentFolder,'MultiGrid'); 

addpath(genpath(pathInput));
addpath(genpath(pathAssembly));
addpath(genpath(pathErrors));
addpath(genpath(pathMeshGeneration));
addpath(genpath(pathFEspace));
addpath(genpath(pathPostProcessing));
addpath(genpath(pathMultiGrid));

%% discontinuos Galerkin - obtain matrices A,b
% to set:
    % problem name (test)
    % domain [0,X] x [0,T]
    % number of finite elements (NT, NX)
    % damping (damp)
    % meshname
    % formulation 0: IP, 1: IPH 
    % dg stability coeff (mu)

test='Test11';   %from InputData %Test11 is the one of the report
Data=DataTest(test);
Data.Degree=2; 
formulation=0; %0: IP, 1: IPH 


plot_sol=1;  %if 1 the solution is plotted 
% (if plot_sol = 1). If there is not an exact solution in Data,
 %   then the analytical plot doesn't make sense

[A,b,femregion,Data] = XT_DG_run(Data,formulation);
[Aj,bj] = blockalize(A,b);

um = zeros(size(Aj,1),1);

for n = 0:10
    um = um + VCycle(Aj,bj,zeros(size(Aj,1),1),5);
    
end