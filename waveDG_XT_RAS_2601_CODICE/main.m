% #######################################################################

% XT-DG WAVE_1D + MULTI_GRID

% #######################################################################

% -------------------MAIN----------------------------- 
% sections:
    % addpath
    % assembly dg matrix
    % multigrid
    % save results

% ----------------------------------------------------

%% addpath
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

% mesh in MeshErrorAnalysis
mesh='ProvaMONO_0105_20100_el.mat';  % domain [0,1]x[0,5], 20 elem in space, 100 in time
Data.X=1;
Data.T=5;
Data.NT =20; 
Data.NX = 100;
Data.damp=0; % DAMPING

formulation=0; %0: IP, 1: IPH 
mu=400;
alpha=mu*(Data.X/Data.NX)/(Data.Degree*Data.Degree); %do not touch 
plot_sol=1;  %if 1 the solution is plotted 
% (if plot_sol = 1). If there is not an exact solution in Data,
 %   then the analytical plot doesn't make sense
Data.meshfile=mesh;
[A,b,femregion,Solutions, ERROR, INFO] = XT_DG_run(Data,alpha,formulation,plot_sol);




%% multigrid

% Generate the block matrixes of A and b
[Aj,bj] = blockalize(A,b);





%% Save results

%--------------------------------------------------------------------------
% --------possibility to save dg matrices in a separate folder ------------

% create name of DG discretized problem
%if formulation ==0
%    f='IP';
%else
%    f='IPH';
%end
%if Data.damp>0
%    name=strcat('mesh_0',num2str(Data.X),'0',num2str(Data.T),'_',num2str(Data.NX),num2str(Data.NT),'_',f,'_mu_',num2str(mu),'_damp');
%else
%    name=strcat('mesh_0',num2str(Data.X),'0',num2str(Data.T),'_',num2str(Data.NX),num2str(Data.NT),'_',f,'_mu_',num2str(mu));
%end
% %% save A and info from the test above
% pathname=fullfile(pwd,'folder_to_create\');

% matfile=fullfile(pathname, name);
% save(matfile,'A','b','Data','femregion','formulation','mesh','Solutions');

% %% load A (if previously saved)
% load(name)
% % back to previous folder
%--------------------------------------------------------------------------

