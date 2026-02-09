function rxm_IdealFlow_CF_Far_Video

% This function generates a video containing the instant flow streamlines
% of the idealized model's near flow.

clc

%% Data pre-processing.
[Xhead_up_BF,Yhead_up_BF,Zhead_up_BF, Xhead_down_BF,Yhead_down_BF,Zhead_down_BF] = Generate_Virtual_Head_UpDownSides;  %[M1 M2]=size(xup1);
xneck = max(Xhead_up_BF(:));
load IdealWF_alpha0.6_Tfine.mat;
Xtail_BF = x(:,1:end-1);
Ytail_BF = y(:,1:end-1);
Ztail_BF = z(:,1:end-1);
Xtail_BF = Xtail_BF+xneck; %size(x)=[ns nt]
clearvars -except Xtail_BF Ytail_BF Ztail_BF Xhead_up_BF Yhead_up_BF Zhead_up_BF  Xhead_down_BF Yhead_down_BF Zhead_down_BF


% Get comoving-frame flow field.
load FarFlow_CF_alpha0.6_Tfine.mat; %
xg=linspace(xl_field(1),xl_field(2),Nx_field);
yg=linspace(yl_field(1),yl_field(2),Ny_field); 
zg=linspace(zl_field(1),zl_field(2),Nz_field); 
[Xg,Yg,Zg]=meshgrid(xg,yg,zg);  
 
% Get comoving-frame sperm shape.
load XNodes_Lab_alpha0.6_Tfine.mat; %
Q=Nhh+Ns;       
%nt_XNodes = size(Xs,2);
Xtail_raw = Xs(Nhh+1:Q,:);  
Ytail_raw = Xs(Q+Nhh+1:2*Q,:); 
Ztail_raw = Xs(2*Q+Nhh+1:3*Q,:);   
% Get the mean flagellar shape direction, based on which the instant sperm shape relative to the comoving frame can be obtained later by rotation.
nt_flow = find(tRange==tRange_flow(end));
Xtail_0 =  Xtail_raw(:,1:nt_flow);
Ytail_0 =  Ytail_raw(:,1:nt_flow);
Ztail_0 =  Ztail_raw(:,1:nt_flow);  
Xtail_1=Xtail_0-repmat(Xtail_0(1,:),Ns,1);
Ytail_1=Ytail_0-repmat(Ytail_0(1,:),Ns,1);
Ztail_1=Ztail_0-repmat(Ztail_0(1,:),Ns,1);   
Xtail_mean= mean(Xtail_1,2); %Ns*1  %The mean flagellar shape relative to the lab frame.
Ytail_mean= mean(Ytail_1,2);
Ztail_mean= mean(Ztail_1,2);
dir_new = [-1, 0, 0];
dir_raw=[Xtail_mean(1)-Xtail_mean(end),Ytail_mean(1)-Ytail_mean(end),Ztail_mean(1)-Ztail_mean(end)];



%% Generate the video. 

% Determine the time period for the video.    
nt_bc = 100;
nt_rev = 600; % This is because the '_Tfine' data have 100 time points for each beat cycle, and the case 'alpha0.6' is of the comoving revolution cycle with 6 beat cycles.

fs = 32;  
fig = 1;
figure(fig) 
v=VideoWriter('FarFlow_alpha0.6_test1','MPEG-4');
v.FrameRate=20;  
open(v);     
for i_nt = 1%:nt_rev+1
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
    
    % Get the instant sperm shape.
    nt0_X = find(tRange==tRange_flow(i_nt));
    % Aligned head.       
    head_tangent_nt0=head_tangent(nt0_X,:)';
    head_normal_nt0=head_normal(nt0_X,:)';
    head_binormal_nt0=cross(head_tangent_nt0,head_normal_nt0);
    B_nt0=[head_tangent_nt0,head_normal_nt0,head_binormal_nt0]; 
    [Bx_nt0_align,By_nt0_align,Bz_nt0_align] = rxm_dir_align(dir_raw,dir_new,B_nt0(1,:)',B_nt0(2,:)',B_nt0(3,:)');
    B_nt0_aligned = [Bx_nt0_align'; By_nt0_align'; Bz_nt0_align']; %3*3.  
    [xup1,xup2,xup3,xdown1,xdown2,xdown3,Ind_xup_neck] = Get_CF_Aligned_Head(B_nt0_aligned);  % xup1: M1*M2.   
    Neck =  [xup1(Ind_xup_neck) xup2(Ind_xup_neck) xup3(Ind_xup_neck)]; 
    % Aligned tail.
    Xtail_nt0 = Xtail_1(:,nt0_X); 
    Ytail_nt0 = Ytail_1(:,nt0_X); 
    Ztail_nt0 = Ztail_1(:,nt0_X); 
    [Xtail_align,Ytail_align,Ztail_align] = rxm_dir_align(dir_raw,dir_new,Xtail_nt0,Ytail_nt0,Ztail_nt0); %Xtail_align: Ns*1.
    Xtail_CF = Xtail_align + Neck(1);  
    Ytail_CF = Ytail_align + Neck(2);
    Ztail_CF = Ztail_align + Neck(3); 
    
    
    % Plot the 3D comoving-frame flow.
    clf;    
    set(gca,'position',[0.09 0.3 0.7 0.7]);  %get(gca,'position')    
    hold on; set(gcf,'color','w'); axis equal; box on; box off; axis off;
    view([-35 20]); % for idealized model's far average flow. 
    XL=[xl_field(1)-2.2 xl_field(2)+2.6]; YL=[yl_field(1)-1.7 yl_field(2)+1.8]; ZL=[zl_field(1)-2.2 zl_field(2)+2.2];
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
    [sx_3D,sy_3D,sz_3D] = meshgrid(xg(1:2:end),yg(21),zg(1:2:end)); 
    Generate_FlowStreamlines(Xg,Yg,Zg,Ug_nt0,Vg_nt0,Wg_nt0,sx_3D,sy_3D,sz_3D, vel_nd);       
    % cross-sectional flow at z=0.
    [sx_3D,sy_3D,sz_3D] = meshgrid(xg(1:2:end),yg(1:2:end),zg(21)); 
    Generate_FlowStreamlines(Xg,Yg,Zg,Ug_nt0,Vg_nt0,Wg_nt0,sx_3D,sy_3D,sz_3D, vel_nd);       
    shading interp
    % Arrows indicating the flow directions.
    if 1 % arrows located at the ends of the flow field. ==> clear visualization.
        [sx_cone1,sy_cone1,sz_cone1] = meshgrid(xg([23]),yg([3 21 40]),zg([2 18 40]));  % arrows located at the top.     
        [sx_cone2,sy_cone2,sz_cone2] = meshgrid(xg([23]),yg([3 21 40]),zg([2 18 40]));  % arrows located at the bottom.     
        [sx_cone3,sy_cone3,sz_cone3] = meshgrid(xg([10 23 40]),yg([3 18 39]),zg([21]));  % arrows located at the left. 
        [sx_cone4,sy_cone4,sz_cone4] = meshgrid(xg([10 23 40]),yg([3 18 39]),zg([21]));  % arrows located at the right.
        [sx_cone5,sy_cone5,sz_cone5] = meshgrid(xg([2 20 43]),yg([21]),zg([8 21 40]));  % arrows located in the front.
        [sx_cone6,sy_cone6,sz_cone6] = meshgrid(xg([2 20 43]),yg([21]),zg([8 21 40]));  % arrows located at the back.       
    else  % arrows located within the flow field. ==> accurate flow directions.
        [sx_cone1,sy_cone1,sz_cone1] = meshgrid(xg([23]),yg([3 21 40]),zg([5 18 37]));  % arrows located at the top.     
        [sx_cone2,sy_cone2,sz_cone2] = meshgrid(xg([23]),yg([3 21 40]),zg([5 18 37]));  % arrows located at the bottom.     
        [sx_cone3,sy_cone3,sz_cone3] = meshgrid(xg([10 23 40]),yg([5 18 37]),zg([21]));  % arrows located at the left. 
        [sx_cone4,sy_cone4,sz_cone4] = meshgrid(xg([10 23 40]),yg([5 18 37]),zg([21]));  % arrows located at the right.
        [sx_cone5,sy_cone5,sz_cone5] = meshgrid(xg([5 20 39]),yg([21]),zg([8 21 40]));  % arrows located in the front.
        [sx_cone6,sy_cone6,sz_cone6] = meshgrid(xg([5 20 39]),yg([21]),zg([8 21 40]));  % arrows located at the back.       
    end
    % arrows located at the top.
    sx_cone1([1 3],:,:)=NaN; sy_cone1([1 3],:,:)=NaN; sz_cone1([1 3],:,:)=NaN;
    sx_cone1(2,:,[1 2])=NaN; sy_cone1(2,:,[1 2])=NaN; sz_cone1(2,:,[1 2])=NaN;
    % arrows located at the bottom.
    sx_cone2([1 3],:,:)=NaN; sy_cone2([1 3],:,:)=NaN; sz_cone2([1 3],:,:)=NaN;
    sx_cone2(2,:,[2 3])=NaN; sy_cone2(2,:,[2 3])=NaN; sz_cone2(2,:,[2 3])=NaN;
    % arrows located at the left.
    sx_cone3(:,[1 3],:)=NaN; sy_cone3(:,[1 3],:)=NaN; sz_cone3(:,[1 3],:)=NaN;
    sx_cone3([2 3],2,:)=NaN; sy_cone3([2 3],2,:)=NaN; sz_cone3([2 3],2,:)=NaN;
    % arrows located at the right.
    sx_cone4(:,[1 3],:)=NaN; sy_cone4(:,[1 3],:)=NaN; sz_cone4(:,[1 3],:)=NaN;
    sx_cone4([1 2],2,:)=NaN; sy_cone4([1 2],2,:)=NaN; sz_cone4([1 2],2,:)=NaN;
    % arrows located in the front.
    sx_cone5(:,:,[1 3])=NaN;  sy_cone5(:,:,[1 3])=NaN;   sz_cone5(:,:,[1 3])=NaN;
    sx_cone5(:,[2 3],2)=NaN;  sy_cone5(:,[2 3],2)=NaN;   sz_cone5(:,[2 3],2)=NaN;
    % arrows located at the back.
    sx_cone6(:,:,[1 3])=NaN;  sy_cone6(:,:,[1 3])=NaN;   sz_cone6(:,:,[1 3])=NaN;
    sx_cone6(:,[1 2],2)=NaN;  sy_cone6(:,[1 2],2)=NaN;   sz_cone6(:,[1 2],2)=NaN;
    % arrow size.    
    Ug_cone = Ug_nt0./vel_raw;  Vg_cone = Vg_nt0./vel_raw;  Wg_cone = Wg_nt0./vel_raw;
    scal = 3; 
    Plot_StreamCone(fig,Xg,Yg,Zg,Ug_cone,Vg_cone,Wg_cone,sx_cone1,sy_cone1,sz_cone1,scal)
    Plot_StreamCone(fig,Xg,Yg,Zg,Ug_cone,Vg_cone,Wg_cone,sx_cone2,sy_cone2,sz_cone2,scal)
    Plot_StreamCone(fig,Xg,Yg,Zg,Ug_cone,Vg_cone,Wg_cone,sx_cone3,sy_cone3,sz_cone3,scal)
    Plot_StreamCone(fig,Xg,Yg,Zg,Ug_cone,Vg_cone,Wg_cone,sx_cone4,sy_cone4,sz_cone4,scal)
    Plot_StreamCone(fig,Xg,Yg,Zg,Ug_cone,Vg_cone,Wg_cone,sx_cone5,sy_cone5,sz_cone5,scal)
    Plot_StreamCone(fig,Xg,Yg,Zg,Ug_cone,Vg_cone,Wg_cone,sx_cone6,sy_cone6,sz_cone6,scal) 
    camlight right
    lighting gouraud
    % Plot the colorbar.
    c=colorbar; 
    c.Position=[0.67 0.35 0.005 0.57];  
    caxis([0 0.0001]);      
    %c.Label.String = 'Flow velocity'; 
    text(15,-17.5,-4,{'Flow velocity'},'FontName','times','FontSize',fs,'rotation',90);  %'FontWeight','bold'  
    c.Ticks=[0 0.0001]; c.TickLabels={'0','$10^{-4}$'}; 
    c.FontName = 'times'; c.TickLabelInterpreter='latex'; c.FontSize=fs; 
    % Plot 2D slice planes.
    if 1 
        % X planes.
        XL1_mat = 0.5;  %[...]
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
        text(-12,11.5,0,{'I'},'FontName','times','interpreter','latex','FontSize',fs);  %xy-plane
        text(-12,0,-12,{'II'},'FontName','times','interpreter','latex','FontSize',fs);  %xz-plane
        text(0,11.5,12,{'III'},'FontName','times','interpreter','latex','FontSize',fs);  %xz-plane  
        
    end    


%% Inset 1: the XYZ coordinate system for the main plot (the 3d flow).
%================================================================================
%================================================================================
    ax1 = axes('Position',[0.08 0.55 0.1 0.1]);  
    view([-35 20]); axis equal; box on; axis off;
    xlabel('x'); ylabel('y'); zlabel('z');
    % Plot Cartesian axes.
    O =[1 0 0];  Cx=[-1; 0; 0];  Cy=[0; -1; 0];  Cz=[0; 0; 1];  C=[Cx Cy Cz];
    Plot_Cartesian_Axes_V1(fig,O,C) 
    camlight right
    lighting gouraud
    % Text: axes labels.
    text(0.4,0,0,{'$x_c$'},'FontName','times','interpreter','latex','FontSize',fs,'FontWeight','bold');   
    text(1,-0.43,0,{'$y_c$'},'FontName','times','interpreter','latex','FontSize',fs,'FontWeight','bold');   
    text(1,0.1,0.52,{'$z_c$'},'FontName','times','interpreter','latex','FontSize',fs,'FontWeight','bold');   
    % Text: varying time moment. 
    t_nd_flow_nt0 = tRange_flow(i_nt)/(2*pi);
    instant_nt0 = roundn(t_nd_flow_nt0,-3);
    text(0.5,0.3,1.3,{['Time: ',num2str(instant_nt0)]},'FontName','times','FontSize',fs+10,'FontWeight','bold');  % time is in units of 'beat cycles'.  
    
    
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
    hold on; set(gcf,'color','w'); axis equal; view([-90 0]); axis off; 
    xlim([xl_field(1)-2 xl_field(2)+2]); ylim([yl_field(1)-2 yl_field(2)+2]); zlim([zl_field(1)-2 zl_field(2)+2]); 
    text(0,0,-12,{'$y_c$'},'FontName','times','interpreter','latex','FontSize',fs);   
    text(0,14,0,{'$z_c$'},'FontName','times','interpreter','latex','FontSize',fs);   
    text(0,17,10,{'III'},'FontName','times','interpreter','latex','FontSize',fs);   
    % Plot streamslice.
    ss=streamslice(Xg,Yg,Zg,Ug_nt0,Vg_nt0,Wg_nt0,[0.5],[],[],0.5);
    set(ss,'color','white','linewidth',3);   
    slice(Xg,Yg,Zg,vel_nd,[0.5],[],[]);
    shading interp; 
    colormap(gca,'parula'); caxis([0 0.0001]);
    % Plot sperm.
    plot3(Xtail_CF,Ytail_CF,Ztail_CF,'color',[36 100 171]./255,'linewidth',2);    %'color',[0.8500    0.3250    0.0980]
    CO_up(:,:,1) = 0.6*ones(size(xup1)); CO_up(:,:,2) = 0*ones(size(xup1)); CO_up(:,:,3) = 0*ones(size(xup1));
    CO_down(:,:,1) = 0*ones(size(xdown1)); CO_down(:,:,2) = 0.6*ones(size(xdown1)); CO_down(:,:,3) = 0*ones(size(xdown1));
    s1= surf(xup1,xup2,xup3,CO_up);
    s2= surf(xdown1,xdown2,xdown3,CO_down);
    s1.EdgeColor='none' ;
    s2.EdgeColor='none' ;
    % Plot XYZ.
    O =[0 10 -10];  Cx=[-21; 0; 0]; Cy=[0; -21; 0]; Cz=[0; 0; 21]; C=[Cx Cy Cz];
    Plot_Cartesian_Axes_V3(fig,O,C)  
    
    % XZ plane projection.
    ax22 = axes('Position',[ps_2d_x ps_2d_y+ps_2d_h+ps_2d_sp ps_2d_w ps_2d_h]);
    hold on; set(gcf,'color','w'); axis equal; view([-180 0]); axis off;
    xlim([xl_field(1)-2 xl_field(2)+2]); ylim([yl_field(1)-2 yl_field(2)+2]); zlim([zl_field(1)-2 zl_field(2)+2]); 
    text(0.5,-1,-12,{'$x_c$'},'FontName','times','interpreter','latex','FontSize',fs);   
    text(15,-1,0,{'$z_c$'},'FontName','times','interpreter','latex','FontSize',fs);   
    text(16,0,10,{'II'},'FontName','times','interpreter','latex','FontSize',fs);  
    % Plot streamslice.
    ss=streamslice(Xg,Yg,Zg,Ug_nt0,Vg_nt0,Wg_nt0,[],[0],[],0.5);
    set(ss,'color','white','linewidth',3);   
    slice(Xg,Yg,Zg,vel_nd,[],[0],[]);
    shading interp; 
    colormap(gca,'parula'); caxis([0 0.0001]);
    % Plot sperm.
    plot3(Xtail_CF,Ytail_CF,Ztail_CF,'color',[36 100 171]./255,'linewidth',2);    %'color',[0.8500    0.3250    0.0980]
    CO_up(:,:,1) = 0.6*ones(size(xup1)); CO_up(:,:,2) = 0*ones(size(xup1)); CO_up(:,:,3) = 0*ones(size(xup1));
    CO_down(:,:,1) = 0*ones(size(xdown1)); CO_down(:,:,2) = 0.6*ones(size(xdown1)); CO_down(:,:,3) = 0*ones(size(xdown1));
    s1= surf(xup1,xup2,xup3,CO_up);
    s2= surf(xdown1,xdown2,xdown3,CO_down);
    s1.EdgeColor='none' ;
    s2.EdgeColor='none' ;
    % Plot XYZ.
    O =[11 1 -10];  Cx=[-22; 0; 0]; Cy=[0; -21; 0]; Cz=[0; 0; 21]; C=[Cx Cy Cz];
    Plot_Cartesian_Axes_V3(fig,O,C)  
    
    % XY plane projection.
    ax23 = axes('Position',[ps_2d_x ps_2d_y+2*(ps_2d_h+ps_2d_sp) ps_2d_w ps_2d_h]);
    hold on; set(gcf,'color','w'); axis equal; view([-180 90]); axis off;
    xlim([xl_field(1)-2 xl_field(2)+2]); ylim([yl_field(1)-2 yl_field(2)+2]); zlim([zl_field(1)-2 zl_field(2)+2]); 
    text(0.5,12,0.5,{'$x_c$'},'FontName','times','interpreter','latex','FontSize',fs);   
    text(15,0,0.5,{'$y_c$'},'FontName','times','interpreter','latex','FontSize',fs);   
    text(15,-10,0,{'I'},'FontName','times','interpreter','latex','FontSize',fs);   
    % Plot streamslice.
    ss=streamslice(Xg,Yg,Zg,Ug_nt0,Vg_nt0,Wg_nt0,[],[],[0],0.5);
    set(ss,'color','white','linewidth',3);   
    slice(Xg,Yg,Zg,vel_nd,[],[],[0]);
    shading interp; 
    colormap(gca,'parula'); caxis([0 0.0001]);
    % Plot sperm.
    plot3(Xtail_CF,Ytail_CF,Ztail_CF,'color',[36 100 171]./255,'linewidth',2);    %'color',[0.8500    0.3250    0.0980]
    CO_up(:,:,1) = 0.6*ones(size(xup1)); CO_up(:,:,2) = 0*ones(size(xup1)); CO_up(:,:,3) = 0*ones(size(xup1));
    CO_down(:,:,1) = 0*ones(size(xdown1)); CO_down(:,:,2) = 0.6*ones(size(xdown1)); CO_down(:,:,3) = 0*ones(size(xdown1));
    s1= surf(xup1,xup2,xup3,CO_up);
    s2= surf(xdown1,xdown2,xdown3,CO_down);
    s1.EdgeColor='none' ;
    s2.EdgeColor='none' ;
    % Plot XYZ.
    O =[11 10 1];  Cx=[-22; 0; 0]; Cy=[0; -21; 0]; Cz=[0; 0; 21]; C=[Cx Cy Cz];
    Plot_Cartesian_Axes_V3(fig,O,C)  
    
    

%% Inset 3: sperm beats in the body and comoving frames, respectively.
%================================================================================
%================================================================================
  
    % Sperm beat in the body frame.
    ax311 = axes('Position',[0.1 0 0.3 0.28]);  
    hold on; set(gcf,'color','w'); axis equal; box on; view(-40,15); axis off;   
    XL = [-0.05,0.95];YL = [-0.18,0.18];ZL = [-0.11,0.11];
    %[min(Xtail_BF(:)) max(Xtail_BF(:)) min(Ytail_BF(:)) max(Ytail_BF(:)) min(Ztail_BF(:)) max(Ztail_BF(:))]=[0.0444    0.9333   -0.1732    0.1732   -0.1054    0.1054]
    xlim([XL(1),XL(2)]);ylim([YL(1),YL(2)]);zlim([ZL(1),ZL(2)]);
    % Plot the tail and head.
    i_nt_WF = mod(i_nt,nt_bc);
    if i_nt_WF==0
        plot3(Xtail_BF(:,end),Ytail_BF(:,end),Ztail_BF(:,end),'color',[36 100 171]./255,'linewidth',2)
    else
        plot3(Xtail_BF(:,i_nt_WF),Ytail_BF(:,i_nt_WF),Ztail_BF(:,i_nt_WF),'color',[36 100 171]./255,'linewidth',2)
    end
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
    ax312 = axes('Position',[0.05 0.11 0.1 0.1]);       
    hold on; set(gcf,'color','w'); axis equal; box on; view(-40,15); axis off; 
    O =[0 0 0]; Cx=[1; 0; 0]; Cy=[0; 1; 0]; Cz=[0; 0; 1]; C=[Cx Cy Cz];
    Plot_Cartesian_Axes_V2(fig,O,C) 
    % Text: axes labels.
    text(1.35,0,0.2,{'$\xi_3$'},'FontName','times','interpreter','latex','FontSize',fs-5);    
    text(-0.3,1.7,0.1,{'$\xi_1$'},'FontName','times','interpreter','latex','FontSize',fs-5);  
    text(-0.1,0,1.65,{'$\xi_2$'},'FontName','times','interpreter','latex','FontSize',fs-5);   
    % Text: sub-title.
    text(0,1.5,2.5,{'Body frame beat'},'FontName','times','FontSize',fs,'FontWeight','bold');   
    camlight right
    lighting gouraud  
    
  
    % Sperm beat in the comoving frame.
    ax321 = axes('Position',[0.47  0.01  0.3  0.28]);   
    hold on; set(gcf,'color','w');axis equal; box on; box off; axis off;
    view([-35 20]);
    XL = [-0.05 0.95];  YL = [-0.14 0.14];  ZL = [-0.14 0.14]; 
    %[min(Xtail_CF(:)) max(Xtail_CF(:)) min(Ytail_CF(:)) max(Ytail_CF(:)) min(Ztail_CF(:)) max(Ztail_CF(:))]=[0.0439    0.9426   -0.1336    0.1360   -0.1373    0.1337] 
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
    ax322 = axes('Position',[0.41  0.05  0.1 0.1]);     
    hold on; set(gcf,'color','w');axis equal; box on; box off; axis off;
    view([-35 20]);
    xlabel('x'); ylabel('y'); zlabel('z');
    O =[0.1 -0.3 -0.25];  Cx=[-1; 0; 0]; Cy=[0; -1; 0]; Cz=[0; 0; 1]; C=[Cx Cy Cz];
    Plot_Cartesian_Axes_V2(fig,O,C) 
    % Text: axes labels.
    text(-1.85,-0.3,-0.25,{'$x_c$'},'FontName','times','interpreter','latex','FontSize',fs,'FontWeight','bold');  
    text(0.1,-1.7,-0.25,{'$y_c$'},'FontName','times','interpreter','latex','FontSize',fs,'FontWeight','bold');   
    text(0,0,1.3,{'$z_c$'},'FontName','times','interpreter','latex','FontSize',fs,'FontWeight','bold');
    % Text: sub-title.
    text(-1.5,0,3.7,{'Comoving frame beat'},'FontName','times','FontSize',fs,'FontWeight','bold');   
    camlight right
    lighting gouraud  
    
 
%% Save the video.
    drawnow;
    frame(i_nt)=getframe(gcf);
    writeVideo(v,frame(i_nt));       
end


