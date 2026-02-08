function rxm_IdealFlowCF_Data
% This function calculates the comoving-frame flow around the idealized sperm, where the sperm head center is located at the origin point.

% We also give PCA decomposition for near-field flow, but not for far-field flow due to its huge computational cost.

clc

if 1 % Calcualte the flow field data.
    
%% Step1: define field size and discretization.

FieldRange = 'n'; % near-filed
%FieldRange = 'f'; % far-filed

if FieldRange == 'n'
    xl_field=[-0.5 1.5]; yl_field=[-1 1]; zl_field=[-1 1];
    Nx_field=21; Ny_field=21; Nz_field=21;   
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



%% Step2: find the 'BCPR', number of beat cycles per trace revolution.

%//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
load XNodes_Lab_alpha0.6_Tfine.mat;
% Trace the distal flagellar point
np = nhh+ns; 
x0 = xs(np,:)';  
y0 = xs(2*np,:)'; 
z0 = xs(3*np,:)'; 
nt = length(x0);
n_beat = tRange/(2*pi); n_beat = roundn(n_beat,-3);
np_beat = find(n_beat==1)-1;


opt = 2; %//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
if opt==1 % alpha=0. 
    np_rev = nt;
elseif opt==2 % 0<alpha<1. 
    np_rev = Fn_TraceRev(x0,y0,z0,np_beat);
elseif opt==3 % alpha=1.
    [np_rev, ~] = fn_rev_TrajType_PerfectH(x0,y0,z0,np_beat);
end

BCPR = np_rev/np_beat % number of beat cycles per trace revolution.


%% Step3: set the position of the flow field relative to the sperm body in the comoving frame.

nt = np_beat*BCPR;
Q=Nhh+Ns;
tRange_flow = tRange(1:nt);


Xtail_0 = Xs(Nhh+1:Q,1:nt); Xtail_1=Xtail_0-repmat(Xtail_0(1,1:nt),Ns,1);
Ytail_0 = Xs(Q+Nhh+1:2*Q,1:nt); Ytail_1=Ytail_0-repmat(Ytail_0(1,1:nt),Ns,1);
Ztail_0 = Xs(2*Q+Nhh+1:3*Q,1:nt); Ztail_1=Ztail_0-repmat(Ztail_0(1,1:nt),Ns,1);

xtail_mean= mean(Xtail_1,2); %Ns*1  
ytail_mean= mean(Ytail_1,2);
ztail_mean= mean(Ztail_1,2);

dir_raw=[-1 0 0];
dir_new=[xtail_mean(1)-xtail_mean(end),ytail_mean(1)-ytail_mean(end),ztail_mean(1)-ztail_mean(end)];
[Xg_align,Yg_align,Zg_align] = rxm_dir_align(dir_raw,dir_new,Xg(:),Yg(:),Zg(:));


xfield_align = [Xg_align(:); Yg_align(:); Zg_align(:)];
Xfield_align = nan(length(xfield_align),nt);
for i_nt=1:nt   
    x0=head_x0(i_nt,:);
    xfield_temp=TranslatePoints(xfield_align,x0);
    Xfield_align(:,i_nt)=xfield_temp;
end




%% Step4: Calculate the flow field.


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
        u=EvaluateVelocityFromForce(Xfield_align(:,i_nt),X,x,f(:,i_nt),epsilon,domain,blockSize);
        [ug,vg,wg]=ExtractComponents(u);

        Ug(:,i_nt) = ug; Vg(:,i_nt) = vg; Wg(:,i_nt) = wg;
    end    
   
    
   if opt==1 % alpha=0. 
       BCPR = NAN;
   end 
   
%//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
save('NearFlow_CF_alpha0.6_Tfine.mat',... 
    'xl_field','yl_field','zl_field','Nx_field','Ny_field','Nz_field',...
    'Ug','Vg','Wg','tRange_flow','BCPR');


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


















