function rxm_ReconFlow_CF_Near_Video_V3

% This function generates a video containing the instant flow streamlines
% of the reconstructed human sperm's near flow.

clc

%% Data pre-processing.

sp=8;
load 01_lab_frame_raw_traces_2017_2018_rxmSSppxyz.mat; 
arc_mat = arclength(X{sp},Y{sp},Z{sp});
arc_mean = mean(arc_mat(end,:));
load FreeSperm_Frequency.mat;  
scale_vel = arc_mean*HF_Freq(sp); % in units o '\mu m/s'.
load NearFlow_CF_sp8_NoBound_Tfine.mat; %
t_dim_flow = t_nd_flow/(2*pi)/HF_Freq(sp);
clearvars -except sp scale_vel t_dim_flow


% Process the body-frame waveform data.
load 03_Head_Fixed_Frame_Method_Rotations_Traces_raw_5_Microns_Difraction_head_spin_rxmSSppxyz.mat;
Xtail = X{sp}; 
Ytail = Y{sp}; 
Ztail = Z{sp}; 
nt=size(Xtail,2);
ns=size(Xtail,1);
for i_nt=1:nt
    arc_mat(:,i_nt) = arclength(Xtail(:,i_nt), Ytail(:,i_nt), Ztail(:,i_nt));
end
arc_mean = mean(arc_mat(end,:)); % measured as 46.35.
Xtail = Xtail/arc_mean;
Ytail = Ytail/arc_mean;
Ztail = Ztail/arc_mean;
% Temporally interpolate the original waveform.
dt=1/90; %sampling frequency of the experimental imaging.
tRange=0:dt:(nt-1)*dt;
t0 = linspace(0,tRange(end),nt);
t1 = t_dim_flow; %t1=linspace(0,tRange(end),5*nt);
s=1:1:ns;
Xtail_BF = interp2(t0,s',Xtail,t1,s'); 
Ytail_BF = interp2(t0,s',Ytail,t1,s'); 
Ztail_BF = interp2(t0,s',Ztail,t1,s'); 
% Generate sperm head: ideal sperm model.
[Xhead_up_BF,Yhead_up_BF,Zhead_up_BF, Xhead_down_BF,Yhead_down_BF,Zhead_down_BF] = Get_BF_IdealSperm_Head(1);  %Get_BF_IdealSperm_Head(arc_mean);  %Xhead_up: M1*M2.
clearvars -except sp scale_vel t_dim_flow Xtail_BF Ytail_BF Ztail_BF Xhead_up_BF Yhead_up_BF Zhead_up_BF  Xhead_down_BF Yhead_down_BF Zhead_down_BF


% Get comoving-frame flow field.
load NearFlow_CF_sp8_NoBound_Tfine.mat; %
xg=linspace(xl_field(1),xl_field(2),Nx_field);
yg=linspace(yl_field(1),yl_field(2),Ny_field); 
zg=linspace(zl_field(1),zl_field(2),Nz_field); 
[Xg,Yg,Zg]=meshgrid(xg,yg,zg);  
 
% Get comoving-frame sperm shape.
load XNodes_sp8_NoBoundStokeslet_Tfine.mat; %
% For the datasets used for figure plot, rather than for generating video, the time points of flow data is coarser than those of sperm mobility.
% Find suitable time moment for the corresponding instant flow.
nt_XNodes=find(t_nd==t_nd_flow(end));
Q=Nhh+Ns;   
Xtail_raw = Xs(Nhh+1:Q,1:nt_XNodes);
Ytail_raw = Xs(Q+Nhh+1:2*Q,1:nt_XNodes); 
Ztail_raw = Xs(2*Q+Nhh+1:3*Q,1:nt_XNodes); 
% Spatial smooth of the flagellar shape, to remove the noise arising from the raw experimental data being reconstructed.
ppx = 0.999; % 'ppx=1' means interpolation without cubic spline.
ppy = 0.999;
ppz = 0.999;  
for i_nt = 1:nt_XNodes
   s_temp = arclength(Xtail_raw(:,i_nt), Ytail_raw(:,i_nt), Ztail_raw(:,i_nt));
   Xtail_sm(:,i_nt) = fnval(csaps(s_temp,Xtail_raw(:,i_nt),ppx),s_temp);
   Ytail_sm(:,i_nt) = fnval(csaps(s_temp,Ytail_raw(:,i_nt),ppy),s_temp);
   Ztail_sm(:,i_nt) = fnval(csaps(s_temp,Ztail_raw(:,i_nt),ppz),s_temp);
end  
% Get the mean flagellar shape direction, based on which the instant sperm shape relative to the comoving frame can be obtained later by rotation.
Xtail_0=Xtail_sm-repmat(Xtail_sm(1,:),Ns,1);
Ytail_0=Ytail_sm-repmat(Ytail_sm(1,:),Ns,1);
Ztail_0=Ztail_sm-repmat(Ztail_sm(1,:),Ns,1);   
Xtail_mean= mean(Xtail_0,2); %Ns*1  %The mean flagellar shape relative to the lab frame.
Ytail_mean= mean(Ytail_0,2);
Ztail_mean= mean(Ztail_0,2);
dir_new = [1, 0, 0];
dir_raw=[Xtail_mean(1)-Xtail_mean(end),Ytail_mean(1)-Ytail_mean(end),Ztail_mean(1)-Ztail_mean(end)];
    



%% Generate the video. 

% Determine the time period for the video.    
[~,nt_rev_temp] = min(abs(t_nd_flow-2*pi*Freq_WF/Freq_rev));
nt_rev = nt_rev_temp  %nt_rev=3112.

fs = 32;  %fs = 46;
fig = 1;
figure(fig) 
v=VideoWriter('NearFlow_sp8_NoBound_test2','MPEG-4');
v.FrameRate=50;
open(v);     

for i_nt = 1%:nt_rev
    i_nt

% Main figure: 3D comoving-frame flow.
%================================================================================
%================================================================================
    %Get instant flow.
    Ug_nt0 = reshape(Ug(:,i_nt),size(Xg,1),size(Xg,2),size(Xg,3));
    Vg_nt0 = reshape(Vg(:,i_nt),size(Xg,1),size(Xg,2),size(Xg,3));
    Wg_nt0 = reshape(Wg(:,i_nt),size(Xg,1),size(Xg,2),size(Xg,3));
    vel_raw=(Ug_nt0.*Ug_nt0+Vg_nt0.*Vg_nt0+Wg_nt0.*Wg_nt0).^0.5;
    vel_nd=vel_raw.*2*pi; % Transform the velocity units into 'flagellum length/ beat cycle'.
    vel_dim = vel_nd.*scale_vel;  % in units of '\mu m/s'.
    
    % Get the instant sperm shape.
    nt0_XNodes=find(t_nd==t_nd_flow(i_nt));
    Xtail_nt0 = Xtail_0(:,nt0_XNodes); 
    Ytail_nt0 = Ytail_0(:,nt0_XNodes); 
    Ztail_nt0 = Ztail_0(:,nt0_XNodes); 
    [Xtail_align,Ytail_align,Ztail_align] = rxm_dir_align(dir_raw,dir_new,Xtail_nt0,Ytail_nt0,Ztail_nt0); %Xtail_align: Ns*1.
    head_tangent_nt0=head_tangent(nt0_XNodes,:)';
    head_normal_nt0=head_normal(nt0_XNodes,:)';
    head_binormal_nt0=cross(head_tangent_nt0,head_normal_nt0);
    B_nt0=[head_tangent_nt0,head_normal_nt0,head_binormal_nt0];
    [Bx_nt0_align,By_nt0_align,Bz_nt0_align] = rxm_dir_align(dir_raw,dir_new,B_nt0(1,:)',B_nt0(2,:)',B_nt0(3,:)');
    B_nt0_aligned = [Bx_nt0_align'; By_nt0_align'; Bz_nt0_align']; %3*3.  
    [xup1,xup2,xup3,xdown1,xdown2,xdown3,Ind_xup_neck] = Get_CF_Aligned_Head(B_nt0_aligned);  % xup1: M1*M2.   
    Neck =  [xup1(Ind_xup_neck) xup2(Ind_xup_neck) xup3(Ind_xup_neck)]; 
    Xtail_CF = Xtail_align + Neck(1);  
    Ytail_CF = Ytail_align + Neck(2);
    Ztail_CF = Ztail_align + Neck(3);
    
    % Plot the 3D comoving-frame flow.
    clf;    
    set(gca,'position',[0.06 0.3 0.7 0.7]);  %get(gca,'position')   
    hold on; set(gcf,'color','w'); axis equal; box on; box off; axis off;
    view([155 15]);  
    XL=[xl_field(1)-0.2 xl_field(2)+0.2]; YL=[yl_field(1)-0.2 yl_field(2)+0.2]; ZL=[zl_field(1)-0.2 zl_field(2)+0.2];
    xlim([XL(1) XL(2)]); ylim([YL(1) YL(2)]); zlim([ZL(1) ZL(2)]);  
    xlabel('x'); ylabel('y'); zlabel('z');
    % Plot sperm.
    plot3(Xtail_CF,Ytail_CF,Ztail_CF,'color',[36 100 171]./255,'linewidth',2);    %'color',[0.8500    0.3250    0.0980]
    CO_up(:,:,1) = 0.6*ones(size(xup1)); CO_up(:,:,2) = 0*ones(size(xup1)); CO_up(:,:,3) = 0*ones(size(xup1));
    CO_down(:,:,1) = 0*ones(size(xdown1)); CO_down(:,:,2) = 0.6*ones(size(xdown1)); CO_down(:,:,3) = 0*ones(size(xdown1));
    s1= surf(xup1,xup2,xup3,CO_up);
    s2= surf(xdown1,xdown2,xdown3,CO_down);
    s1.EdgeColor='none' ;
    s2.EdgeColor='none' ;
    % Plot flow streamlines.
    % cross-sectional flow at y=0.
    [sx_3D,sy_3D,sz_3D] = meshgrid(xg(:),yg(11),zg(:)); 
    Generate_FlowStreamlines(Xg,Yg,Zg,Ug_nt0,Vg_nt0,Wg_nt0,sx_3D,sy_3D,sz_3D, vel_dim);       
    % cross-sectional flow at z=0.
    [sx_3D,sy_3D,sz_3D] = meshgrid(xg(:),yg(:),zg(11)); 
    Generate_FlowStreamlines(Xg,Yg,Zg,Ug_nt0,Vg_nt0,Wg_nt0,sx_3D,sy_3D,sz_3D, vel_dim);       
    shading interp
    camlight right
    lighting gouraud
    % Plot the colorbar.
    c=colorbar; 
    c.Position=[0.66 0.35 0.005 0.6];  
    caxis([0 5]);      
    c.Label.String = 'Flow velocity (\mum/s)'; 
    c.Ticks=[0 5]; c.TickLabels={'0','5'}; 
    c.FontName = 'times'; c.TickLabelInterpreter='latex'; c.FontSize=fs; 
    % Plot 2D slice planes.
    if 1 
        % X planes.
        XL1_mat = -0.5; %[-4 -1.5 0 1 4];
        for i_x = 1:length(XL1_mat)
            XL_temp=[XL1_mat(i_x) XL(2)];
            YL_temp=[YL(1) YL(2)];
            ZL_temp=[ZL(1) ZL(2)]; 
            [xgrid,ygrid,zgrid] = meshgrid(linspace(XL_temp(1),XL_temp(2),10),...
                    linspace(YL_temp(1),YL_temp(2),10),...
                    linspace(ZL_temp(1),ZL_temp(2),10));     
            CO(:,:,1) = ones(10)*.0; % red 
            CO(:,:,2) = ones(10)*.4470; % green 
            CO(:,:,3) = ones(10)*0.7410; % blue 
            surfX = surf(reshape(xgrid(:,1,:),[10 10]),reshape(ygrid(:,1,:),[10 10]),reshape(zgrid(:,1,:),[10,10]), CO,'FaceAlpha',0.1);  
            surfX.EdgeColor='none';
        end
        % Y plane.
        XL_temp=[XL(1) XL(2)];
        YL_temp=[0 YL(2)];
        ZL_temp=[ZL(1) ZL(2)]; 
        [xgrid,ygrid,zgrid] = meshgrid(linspace(XL_temp(1),XL_temp(2),10),...
                    linspace(YL_temp(1),YL_temp(2),10),...
                    linspace(ZL_temp(1),ZL_temp(2),10));     
        CO(:,:,1) = ones(10)*.0; % red 
        CO(:,:,2) = ones(10)*.4470; % green 
        CO(:,:,3) = ones(10)*0.7410; % blue 
        surfY = surf(reshape(xgrid(1,:,:),[10 10]),reshape(ygrid(1,:,:),[10 10]),reshape(zgrid(1,:,:),[10,10]), CO,'FaceAlpha',0.1);  
        surfY.EdgeColor='none';
        % Z plane.          
        XL_temp=[XL(1) XL(2)];
        YL_temp=[YL(1) YL(2)];
        ZL_temp=[0 ZL(2)]; 
        [xgrid,ygrid,zgrid] = meshgrid(linspace(XL_temp(1),XL_temp(2),10),...
                    linspace(YL_temp(1),YL_temp(2),10),...
                    linspace(ZL_temp(1),ZL_temp(2),10));     
        CO(:,:,1) = ones(10)*.0; % red 
        CO(:,:,2) = ones(10)*.4470; % green 
        CO(:,:,3) = ones(10)*0.7410; % blue 
        surfZ = surf(reshape(xgrid(:,:,1),[10 10]),reshape(ygrid(:,:,1),[10 10]),reshape(zgrid(:,:,1),[10,10]), CO,'FaceAlpha',0.1);  
        surfZ.EdgeColor='none';
        %Texts indicating the three cross sections.
        text(0.65,-1.15,0.03,{'I'},'FontName','times','interpreter','latex','FontSize',fs);  %xy-plane   
        text(0.65,0,1.18,{'II'},'FontName','times','interpreter','latex','FontSize',fs);  %xz-plane  
        text(-0.5,-1.15,1.15,{'III'},'FontName','times','interpreter','latex','FontSize',fs);  %xz-plane   
    end    

    
%% Inset 1: the XYZ coordinate system for the main plot (the 3d flow).
%================================================================================
%================================================================================
    ax1 = axes('Position',[0.08 0.55 0.1 0.1]);  
    view([155 15]); axis equal; box on; axis off;
    xlabel('x'); ylabel('y'); zlabel('z');
    % Plot Cartesian axes.
    O =[1 0 0];  Cx=[1; 0; 0];  Cy=[0; 1; 0];  Cz=[0; 0; 1];  C=[Cx Cy Cz];
    Plot_Cartesian_Axes(fig,O,C) 
    camlight right
    lighting gouraud
    % Text: axes labels.
    text(1.6,0,0,{'$x_c$'},'FontName','times','interpreter','latex','FontSize',fs,'FontWeight','bold');   
    text(1,0.45,0,{'$y_c$'},'FontName','times','interpreter','latex','FontSize',fs,'FontWeight','bold');   
    text(1,-0.1,0.48,{'$z_c$'},'FontName','times','interpreter','latex','FontSize',fs,'FontWeight','bold');   
    % Text: varying time moment. 
    t_dim_flow_nt0 = t_nd_flow(i_nt)/(2*pi)/Freq_WF;
    instant_nt0 = roundn(t_dim_flow_nt0,-3);
    text(1.8,-0.3,1.3,{['Time: ',num2str(instant_nt0),'s']},'FontName','times','FontSize',fs+10,'FontWeight','bold');   
   
    
    
%% Inset 2: the three 2D stream slices.
%================================================================================
%================================================================================
    ps_2d_x = 0.75;  
    ps_2d_y = 0.08;  
    ps_2d_w = 0.23; 
    ps_2d_h = 0.23;  
    ps_2d_sp = 0.085; 
    
    % YZ plane projection.
    ax21 = axes('Position',[ps_2d_x ps_2d_y ps_2d_w ps_2d_h]);
    hold on; set(gcf,'color','w'); axis equal; view([90 0]);  
    xlim([xl_field(1) xl_field(2)]); ylim([yl_field(1) yl_field(2)]); zlim([zl_field(1) zl_field(2)]); 
    xticks([]); yticks([]); zticks([]);
    ylh=ylabel('$y_c$','interpreter','latex','FontSize',fs);
    zlh=zlabel('$z_c$','interpreter','latex','FontSize',fs,'Rotation',0);
    set(ylh,'position',[0,0,-1.05]);  set(zlh,'position',[0,-1.2,0])
    text(0,-1.5,1,{'III'},'FontName','times','interpreter','latex','FontSize',fs);   
    % Plot streamslice.
    ss=streamslice(Xg,Yg,Zg,Ug_nt0,Vg_nt0,Wg_nt0,[-0.5],[],[],0.5);
    set(ss,'color','white','linewidth',3);   
    slice(Xg,Yg,Zg,vel_dim,[-0.5],[],[]);
    shading interp; 
    colormap(gca,'parula'); caxis([0 5]);
    % Plot sperm.
    plot3(Xtail_CF,Ytail_CF,Ztail_CF,'color',[36 100 171]./255,'linewidth',2);    %'color',[0.8500    0.3250    0.0980]
    CO_up(:,:,1) = 0.6*ones(size(xup1)); CO_up(:,:,2) = 0*ones(size(xup1)); CO_up(:,:,3) = 0*ones(size(xup1));
    CO_down(:,:,1) = 0*ones(size(xdown1)); CO_down(:,:,2) = 0.6*ones(size(xdown1)); CO_down(:,:,3) = 0*ones(size(xdown1));
    s1= surf(xup1,xup2,xup3,CO_up);
    s2= surf(xdown1,xdown2,xdown3,CO_down);
    s1.EdgeColor='none' ;
    s2.EdgeColor='none' ;
    
    % XZ plane projection.
    ax22 = axes('Position',[ps_2d_x ps_2d_y+ps_2d_h+ps_2d_sp ps_2d_w ps_2d_h]);
    hold on; set(gcf,'color','w'); axis equal; view([0 0]);  
    xlim([xl_field(1) xl_field(2)]); ylim([yl_field(1) yl_field(2)]); zlim([zl_field(1) zl_field(2)]);  
    xticks([]); yticks([]); zticks([]);
    xlh=xlabel('$x_c$','interpreter','latex','FontSize',fs);
    zlh=zlabel('$z_c$','interpreter','latex','FontSize',fs,'Rotation',0);
    set(xlh,'position',[-0.5,-1,-1.05]);  set(zlh,'position',[-1.7,-1,0])
    text(-1.85,0,1,{'II'},'FontName','times','interpreter','latex','FontSize',fs);  
    % Plot streamslice.
    ss=streamslice(Xg,Yg,Zg,Ug_nt0,Vg_nt0,Wg_nt0,[],[0],[],0.5);
    set(ss,'color','white','linewidth',3);   
    slice(Xg,Yg,Zg,vel_dim,[],[0],[]);
    shading interp; 
    colormap(gca,'parula'); caxis([0 5]);
    % Plot sperm.
    plot3(Xtail_CF,Ytail_CF,Ztail_CF,'color',[36 100 171]./255,'linewidth',2);    %'color',[0.8500    0.3250    0.0980]
    CO_up(:,:,1) = 0.6*ones(size(xup1)); CO_up(:,:,2) = 0*ones(size(xup1)); CO_up(:,:,3) = 0*ones(size(xup1));
    CO_down(:,:,1) = 0*ones(size(xdown1)); CO_down(:,:,2) = 0.6*ones(size(xdown1)); CO_down(:,:,3) = 0*ones(size(xdown1));
    s1= surf(xup1,xup2,xup3,CO_up);
    s2= surf(xdown1,xdown2,xdown3,CO_down);
    s1.EdgeColor='none' ;
    s2.EdgeColor='none' ;
    
    % XY plane projection.
    ax23 = axes('Position',[ps_2d_x ps_2d_y+2*(ps_2d_h+ps_2d_sp) ps_2d_w ps_2d_h]);
    hold on; set(gcf,'color','w'); axis equal; view([0 90]); 
    xticks([]); yticks([]); zticks([]);
    xlim([xl_field(1) xl_field(2)]); ylim([yl_field(1) yl_field(2)]); zlim([zl_field(1) zl_field(2)]);  
    xlh=xlabel('$x_c$','interpreter','latex','FontSize',fs);
    ylh=ylabel('$y_c$','interpreter','latex','FontSize',fs);
    set(xlh,'position',[-0.5,-1.05,0.5]);  set(ylh,'position',[-1.55,0,0.5])
    text(-1.73,1,0,{'I'},'FontName','times','interpreter','latex','FontSize',fs);   
    % Plot streamslice.
    ss=streamslice(Xg,Yg,Zg,Ug_nt0,Vg_nt0,Wg_nt0,[],[],[0],0.5);
    set(ss,'color','white','linewidth',3);   
    slice(Xg,Yg,Zg,vel_dim,[],[],[0]);
    shading interp; 
    colormap(gca,'parula'); caxis([0 5]);
    % Plot sperm.
    plot3(Xtail_CF,Ytail_CF,Ztail_CF,'color',[36 100 171]./255,'linewidth',2);    %'color',[0.8500    0.3250    0.0980]
    CO_up(:,:,1) = 0.6*ones(size(xup1)); CO_up(:,:,2) = 0*ones(size(xup1)); CO_up(:,:,3) = 0*ones(size(xup1));
    CO_down(:,:,1) = 0*ones(size(xdown1)); CO_down(:,:,2) = 0.6*ones(size(xdown1)); CO_down(:,:,3) = 0*ones(size(xdown1));
    s1= surf(xup1,xup2,xup3,CO_up);
    s2= surf(xdown1,xdown2,xdown3,CO_down);
    s1.EdgeColor='none' ;
    s2.EdgeColor='none' ;
    


%% Inset 3: sperm beats in the body and comoving frames, respectively.
%================================================================================
%================================================================================
    
    % Sperm beat in the body frame.
    ax311 = axes('Position',[0.05 0 0.38 0.38]);  
    hold on; set(gcf,'color','w'); axis equal; box on; view(-45,16); axis off;   
    XL = [-0.15 0.95];  YL = [-0.6 0.74];  ZL = [-0.72 0.54];  
    xlim([XL(1),XL(2)]);ylim([YL(1),YL(2)]);zlim([ZL(1),ZL(2)]);
    % Plot the tail and head.
    plot3(Xtail_BF(:,i_nt),Ytail_BF(:,i_nt),Ztail_BF(:,i_nt),'color',[36 100 171]./255,'linewidth',2)
    source_ligth = [90 80];
    k = [0.9 1 1 5]; 
    s1= surfl(Xhead_up_BF,Yhead_up_BF,Zhead_up_BF,source_ligth,k );
    s2= surfl(Xhead_down_BF,Yhead_down_BF,Zhead_down_BF,source_ligth,k );
    s1.EdgeColor='none' ;
    s1.FaceColor=  [0.6 0 0 ];
    s2.EdgeColor='none' ;
    s2.FaceColor=  [0 0.6 0]; 
    camlight right
    lighting gouraud  
    
    
    % Plot the Cartesian axes for the body frame.
    ax312 = axes('Position',[0.09 0.08 0.1 0.1]);   
    hold on; set(gcf,'color','w'); axis equal; box on; view(-45,16);  axis off; 
    O =[0 0 0]; Cx=[1; 0; 0]; Cy=[0; 1; 0]; Cz=[0; 0; 1]; C=[Cx Cy Cz];
    Plot_Cartesian_Axes_V2(fig,O,C) 
    % Text: axes labels.
    text(1.4,0,0,{'$\xi_3$'},'FontName','times','interpreter','latex','FontSize',fs-5);   
    text(-0.3,1.7,0,{'$\xi_1$'},'FontName','times','interpreter','latex','FontSize',fs-5);  
    text(-0.1,0,1.65,{'$\xi_2$'},'FontName','times','interpreter','latex','FontSize',fs-5);   
    % Text: sub-title.
    text(-3,1.5,4,{'Body frame beat'},'FontName','times','FontSize',fs,'FontWeight','bold');  
    camlight right
    lighting gouraud  
    
    
    % Sperm beat in the comoving frame.
    ax321 = axes('Position',[0.42  0.01  0.4  0.28]);   
    hold on; set(gcf,'color','w');axis equal; box on; box off; axis off;
    view([155 15]);
    XL = [-1.1 0.05];  YL = [-0.3 0.3];  ZL = [-0.3 0.3];  
    xlim([XL(1) XL(2)]); ylim([YL(1) YL(2)]); zlim([ZL(1) ZL(2)]);  
    xlabel('x'); ylabel('y'); zlabel('z');
    % Plot sperm.
    plot3(Xtail_CF,Ytail_CF,Ztail_CF,'color',[36 100 171]./255,'linewidth',2);    %'color',[0.8500    0.3250    0.0980]
    CO_up(:,:,1) = 0.6*ones(size(xup1)); CO_up(:,:,2) = 0*ones(size(xup1)); CO_up(:,:,3) = 0*ones(size(xup1));
    CO_down(:,:,1) = 0*ones(size(xdown1)); CO_down(:,:,2) = 0.6*ones(size(xdown1)); CO_down(:,:,3) = 0*ones(size(xdown1));
    s1= surf(xup1,xup2,xup3,CO_up);
    s2= surf(xdown1,xdown2,xdown3,CO_down);
    s1.EdgeColor='none' ;
    s2.EdgeColor='none' ;
    camlight right
    lighting gouraud  
    
    
    % Plot the Cartesian axes for the comoving frame.
    ax322 = axes('Position',[0.43  0.05  0.1 0.1]);    
    hold on; set(gcf,'color','w');axis equal; box on; box off; axis off;
    view([155 15]);
    xlabel('x'); ylabel('y'); zlabel('z');
    O =[0.1 -0.3 -0.25];  Cx=[1; 0; 0]; Cy=[0; 1; 0]; Cz=[0; 0; 1]; C=[Cx Cy Cz];
    Plot_Cartesian_Axes_V2(fig,O,C) 
    % Text: axes labels.
    text(2,-0.3,-0.25,{'$x_c$'},'FontName','times','interpreter','latex','FontSize',fs,'FontWeight','bold');  
    text(0.1,1.1,-0.25,{'$y_c$'},'FontName','times','interpreter','latex','FontSize',fs,'FontWeight','bold');   
    text(0.3,-0.4,1.35,{'$z_c$'},'FontName','times','interpreter','latex','FontSize',fs,'FontWeight','bold');
    % Text: sub-title.
    text(2,0,3.6,{'Comoving frame beat'},'FontName','times','FontSize',fs,'FontWeight','bold');   
    camlight right
    lighting gouraud  
    
    
%% Save the video.
    drawnow;
    frame(i_nt)=getframe(gcf);
    writeVideo(v,frame(i_nt));       
end


