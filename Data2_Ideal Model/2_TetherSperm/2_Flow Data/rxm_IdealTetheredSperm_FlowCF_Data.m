function rxm_IdealTetheredSperm_FlowCF_Data
% This function calculates the comoving-frame flow around the idealized sperm, where the sperm head center is located at the origin point.

% We also give PCA decomposition for near-field flow, but not for far-field flow due to its huge computational cost.

clc

if 1 % Calcualte the flow field data.
    
%% Step1: define field size and discretization.

FieldRange = 'n'; % near-filed
%FieldRange = 'm'; % middle-filed
%FieldRange = 'f'; % far-filed

if FieldRange == 'n'
    xl_field=[-0.5 1.5]; yl_field=[-1 1]; zl_field=[-1 1];
    Nx_field=21; Ny_field=21; Nz_field=21;   
elseif FieldRange == 'm'
    % For 'No Bound' case
    xl_field=[-1 2]; yl_field=[-1.5 1.5]; zl_field=[-0.2 0.2];
    Nx_field=31; Ny_field=31; Nz_field=5;
elseif FieldRange == 'f'
    % For 'No Bound' case
    xl_field=[-10 11]; yl_field=[-10 10]; zl_field=[-10 10];
    Nx_field=43; Ny_field=41; Nz_field=41;
    % For '1 Bound' case
    % bound_h=0.2;
    % xl_field=[-4 5]; yl_field=[-4 4]; zl_field=[-bound_h 1-bound_h];
    % Nx_field=46; Ny_field=41; Nz_field=21;
end

xg=linspace(xl_field(1),xl_field(2),Nx_field);
yg=linspace(yl_field(1),yl_field(2),Ny_field); 
zg=linspace(zl_field(1),zl_field(2),Nz_field); 
% plotting grid to form the local flow field
[Xg,Yg,Zg]=meshgrid(xg,yg,zg);
Xfield = [Xg(:);Yg(:);Zg(:)];


%% Step2: Calculate the flow field.


%load XNodes_TetheredSperm_alpha0.4.mat;
load XNodes_TetheredSperm_alpha0.4_V2.mat;
nt = size(xs,2);

% numerical parameters
    epsilon=0.25/45;
    domain='i';
    blockSize=0.2;

    Ug = nan(length(Xg(:)),nt); 
    Vg = nan(length(Yg(:)),nt); 
    Wg = nan(length(Zg(:)),nt);
    for i_nt=1:nt
        i_nt
        % Assemble allocation points: sperm and boundary
        x=MergeVectorGrids(xs(:,i_nt),xb); 
        X=MergeVectorGrids(Xs(:,i_nt),Xb); 
        
        % evaluate velocity field on grid and then save the data
        u=EvaluateVelocityFromForce(Xfield,X,x,f(:,i_nt),epsilon,domain,blockSize);
        [ug,vg,wg]=ExtractComponents(u);

        Ug(:,i_nt) = ug; Vg(:,i_nt) = vg; Wg(:,i_nt) = wg;
    end    
   

%//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
save('NearFlow_Tethered_alpha0.4_V2.mat',... 
    'xl_field','yl_field','zl_field','Nx_field','Ny_field','Nz_field',...
    'Ug','Vg','Wg','tRange');


end


%% Flow PCA calculation: only for near field.

%//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
% Take the idealized model using alpha=0.6 as the example for PCA analysis.
if 0
    
    load NearFlow_CF_alpha0.6_Tfine.mat;
    np = size(Ug,1);
    %Ug = Ug(:,1:600); Vg = Vg(:,1:600); Wg = Wg(:,1:600);
    [PCA_mode,Time_coef, nmode,cumulative_variance]=Extract_flow_PCA(Ug,Vg,Wg); 
    size(PCA_mode)
    PCA_mode_Ug = PCA_mode(1:np,:); 
    PCA_mode_Vg = PCA_mode(np+1:2*np,:); 
    PCA_mode_Wg = PCA_mode(2*np+1:3*np,:);
    
    save(['NearFlow_CF_alpha0.6_Tfine_PCA_Accuracy99.99.mat'],'PCA_mode_Ug','PCA_mode_Vg','PCA_mode_Wg','Time_coef','nmode','cumulative_variance');  
    %save(['NearFlow_CF_alpha0.6_TfineTshort_PCA_Accuracy99.99.mat'],'PCA_mode_Ug','PCA_mode_Vg','PCA_mode_Wg','Time_coef','nmode','cumulative_variance');  


end


















