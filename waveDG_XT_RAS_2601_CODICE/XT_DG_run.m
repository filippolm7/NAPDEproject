% DG discretization
% IP and IPh formulation

function [A,b,femregion,Data] = XT_DG_run(Data,formulation)


% mesh in MeshErrorAnalysis
mesh='ProvaMONO_0105_20100_el.mat';  
Data.X=1;
Data.T=5;       % domain [0,1]x[0,5]
Data.NT =20;     % Default 20 elem in space, 100 in time
Data.NX = 100;
Data.damp=0; % DAMPING
Data.nqn = 2*Data.Degree + 1;
Data.meshfile=mesh;

mu=400;
alpha=mu*(Data.X/Data.NX)/(Data.Degree*Data.Degree); %do not touch 

Data.fem = Data.Degree;
Data.alp=alpha;
% set up quadrature nodes
Data.nqn = 2*Data.fem + 1;


%% load meshfile
load(Data.meshfile);
disp('Uploading mesh...');
disp('------------------------------------------------------')


Data.domain = Dati.domain;
region.id = E_tag;

%% checking tags 
% keep tag A that comes from older implentation

for i = 1 : length(E_tag)   
    
    for j = 1 : length(Data.tag_ac)
        if E_tag(i) == Data.tag_ac(j)
            region.tag(i,1) = 'A';
        end
    end

end


%% computing neighbourung elements
[neighbour] = neighbours_new(region,Data);


%% plot mesh
plot_poly_mesh(Data,region,neighbour)
hold on;
for i = 1 : size(region.coord,1)
    text(region.coord(i,1),region.coord(i,2), num2str(i));
end


%%  creation of the finite element space
[femregion] = create_dof_new(Data,region);

%%  MATRIX ASSEMBLY

disp('Matrices computation ... ');

if formulation==1
    [MatricesAc]  = matrix2D_iph(femregion,neighbour,Data);
    disp('IPH formulation')
elseif formulation==0
    [MatricesAc]  = matrix2D_ip(femregion,neighbour,Data);
    disp('IP formulation')
else
    disp('formulation error')
end


disp('Done')
disp('------------------------------------------------------')


%% RIGHT HAND SIDE

disp('Computing rhs ... ');
if formulation==1
    [fx_a] = evaluate_f_acu_iph(neighbour,femregion,Data);
elseif formulation==0
    [fx_a] = evaluate_f_acu_ip(neighbour,femregion,Data);
else
    disp('formulation error')
end
disp('Done')
disp('------------------------------------------------------')


%% return
A=MatricesAc.block;
b=fx_a;

end
