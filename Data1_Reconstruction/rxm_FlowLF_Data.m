function rxm_FlowLF_Data

% This function calculates the flow field around a sperm, relative to the
% stationary lab frame of reference.

% The primary purpose of the lab-frame flow (data) is to check our
% correct implementation of the image system by confirming that the
% flow velocity on boundary surface equal to zero. Therefore, we don't need
% full time-evolution of calculation.


clc


%% Step1: define the flow field size and discretization (near field).

xl_field=[-1.5 0.5]; yl_field=[-1 1]; zl_field=[-1 1];
Nx_field=21; Ny_field=21; Nz_field=21;

% 3D grid space of the local flow field.
xg=linspace(xl_field(1),xl_field(2),Nx_field);
yg=linspace(yl_field(1),yl_field(2),Ny_field); 
zg=linspace(zl_field(1),zl_field(2),Nz_field); 
[Xg,Yg,Zg]=meshgrid(xg,yg,zg);

Xfield = [Xg(:); Yg(:); Zg(:)];



%% Step2: Calculate instantaneous flow field (with/ without a boundary).

    % numerical parameters
    epsilon=0.25/45;
    domain='h';    %Blakelet 
    blockSize=0.2;


    load XNodes_sp8_BoundBlakelet_H0.2.mat; 
    nBC_raw = t_nd(end)/(2*pi); %t_nd: non-dimensioanl time.
    nt_raw = length(t_nd);
    dt_beat_raw = nt_raw/nBC_raw;
    dt_beat_temp=20;
    ddt = round(dt_beat_raw/dt_beat_temp);
    nt_flow = 40; % Specify the time point number to be calculated.
    t_nd_flow = t_nd(1:ddt:1+(nt_flow-1)*ddt);
    

    Ug = nan(length(Xg(:)),nt_flow); 
    Vg = nan(length(Yg(:)),nt_flow); 
    Wg = nan(length(Zg(:)),nt_flow);
    
    
    for i_nt=1:nt_flow 
        i_nt
        % Assemble allocation points: sperm and boundary
        x=xs(:,1+(i_nt-1)*ddt);%   x=MergeVectorGrids(xs(:,1+(i_nt-1)*dnp),xb); <-- 'xb=[]'.
        X=Xs(:,1+(i_nt-1)*ddt);%   X=MergeVectorGrids(Xs(:,1+(i_nt-1)*dnp),Xb); <-- 'xb=[]'.
        
        % evaluate velocity field on grid and then save the data
        u=EvaluateVelocityFromForce(Xfield,X,x,f(:,1+(i_nt-1)*ddt),epsilon,domain,blockSize);
        [ug,vg,wg]=ExtractComponents(u);

        Ug(:,i_nt) = ug; 
        Vg(:,i_nt) = vg; 
        Wg(:,i_nt) = wg;
    end
    
    
    
    
save('NearFlow_LF_sp8_BoundH0.2.mat',...
     'xl_field','yl_field','zl_field','Nx_field','Ny_field','Nz_field',...
     'Ug','Vg','Wg','t_nd_flow');





