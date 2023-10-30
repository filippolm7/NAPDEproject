    
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

mesh='ProvaMONO_0105_20100_el.mat';
Data.X=1;
Data.T=5;       % domain [0,1]x[0,5]
Data.NT =20;     % Default 20 elem in space, 100 in time
Data.NX = 100;
Data.damp=0; % DAMPING
Data.nqn = 2*Data.Degree + 1;
Data.meshfile=mesh;

[A,b,femregion,Data] = XT_DG_run(Data,formulation);

[Aj,bj] = blockalize(A,b);


R = RMatrix(size(Aj,1),1);
Ah = R*Aj*R';
bh = R*bj;

uh_wh=Ah\bh;
[uj,wj] = zoop(uh_wh);
Solutions.phi_h=uj;
Solutions.dot_phi_h=wj;

%% scatter plot of the solution 
%if plot_sol==1

  %  [GUp,GWp,Gphi,GUe] = plot_solution_dis(Data,femregion,Solutions,0.1); %t=0.1
%    scatter_plot(GUp,GWp,Gphi,GUe,0,Data); 
%end

um = zeros(size(Aj,1),1);



%for n = 0:10
 %   um = um + VCycle(Aj,bj,zeros(size(Aj,1),1),5);
    
%end

