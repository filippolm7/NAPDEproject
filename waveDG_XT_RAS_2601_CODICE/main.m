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
formulation=0; %0: IP, 1: IPH 


plot_sol=1;  %if 1 the solution is plotted 
% (if plot_sol = 1). If there is not an exact solution in Data,
 %   then the analytical plot doesn't make sense

%% Assemble 

[A,b,femregion,Data] = XT_DG_run(Data,formulation);

%% Solve
uh_wh=A\b;
m=size(uh_wh)/2;
wh=uh_wh((m+1):end);
uh=uh_wh(1:m);
Solutions.phi_h=uh;
Solutions.dot_phi_h=wh;

%% Compute Errors
% missing part not studied in this project

ERROR=-1;
INFO=-1;



%% scatter plot of the solution 
if plot_sol==1

    [GUp,GWp,Gphi,GUe] = plot_solution_dis(Data,femregion,Solutions,0.1); %t=0.1
    scatter_plot(GUp,GWp,Gphi,GUe,0,Data); 
end


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

