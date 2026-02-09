function rxm_ReconFlow_AverageNear_CrossSection_Plot

% This function plots the average flow field around the reconstructed human sperm, relative to the comoving frame of reference.
% Figures:
% 1. 3D average flow. Parameters can be adjusted: with/without a boundary, full/ cross-section streamlines.
% 2. 2D average flow.
% 3. Colorbar.

clc

sp=8;

load 01_lab_frame_raw_traces_2017_2018_rxmSSppxyz.mat; 
arc_mat = arclength(X{sp},Y{sp},Z{sp});
arc_mean = mean(arc_mat(end,:));

load FreeSperm_Frequency.mat;  
scale_vel = arc_mean*HF_Freq(sp); % in units o '\mu m/s'.

clearvars -except scale_vel


%% Fig. 1: 3D flow.

if 1
    %=====================================================================================================
    % Get flow field.
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////
    %load NearFlow_CF_sp8_BoundH1.mat;
    load NearFlow_CF_sp8_NoBound.mat;
    
    xg=linspace(xl_field(1),xl_field(2),Nx_field);
    yg=linspace(yl_field(1),yl_field(2),Ny_field); 
    zg=linspace(zl_field(1),zl_field(2),Nz_field); 
    [Xg,Yg,Zg]=meshgrid(xg,yg,zg);   
    
    Ug_mean = mean(Ug,2);
    Vg_mean = mean(Vg,2);
    Wg_mean = mean(Wg,2);
    Ug_mean = reshape(Ug_mean,size(Xg,1),size(Xg,2),size(Xg,3));
    Vg_mean = reshape(Vg_mean,size(Xg,1),size(Xg,2),size(Xg,3));
    Wg_mean = reshape(Wg_mean,size(Xg,1),size(Xg,2),size(Xg,3));
    
    vel_raw=(Ug_mean.*Ug_mean+Vg_mean.*Vg_mean+Wg_mean.*Wg_mean).^0.5; 
    vel_nd=vel_raw.*2*pi;  %in units of 'flagellum length/ beat cycle',
    vel_dim = vel_nd.*scale_vel;  % in units of '\mu m/s'.
   
    
    
    %=====================================================================================================
    % Get sperm shape.
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////
    %load XNodes_sp8_BoundBlakelet_H1.mat; 
    load XNodes_sp8_NoBoundStokeslet.mat;
    
    % For the datasets used for figure plot, rather than for generating video, the time points of flow data is coarser than those of sperm mobility.
    % Find suitable time period to average: complete trace revolution periods, as selected for the flow-field time.
    nt_XNodes=find(t_nd==t_nd_flow(end));
    Q=Nhh+Ns;   
    Xtail_0 = Xs(Nhh+1:Q,1:nt_XNodes); 
    Ytail_0 = Xs(Q+Nhh+1:2*Q,1:nt_XNodes); 
    Ztail_0 = Xs(2*Q+Nhh+1:3*Q,1:nt_XNodes); 
    % Spatial smooth of the flagellar shape, to remove the noise arising from the raw experimental data being reconstructed.
    ppx = 0.999; % 'ppx=1' means interpolation without cubic spline.
    ppy = 0.999;
    ppz = 0.999;  
    for i_nt = 1:nt_XNodes
        s_temp = arclength(Xtail_0(:,i_nt), Ytail_0(:,i_nt), Ztail_0(:,i_nt));
        Xtail_0(:,i_nt) = fnval(csaps(s_temp,Xtail_0(:,i_nt),ppx),s_temp);
        Ytail_0(:,i_nt) = fnval(csaps(s_temp,Ytail_0(:,i_nt),ppy),s_temp);
        Ztail_0(:,i_nt) = fnval(csaps(s_temp,Ztail_0(:,i_nt),ppz),s_temp);
    end    
    
    % Get average sperm body shape.
    dir_align = [1, 0, 0];
    [Xtail, Ytail, Ztail] = Get_Aligned_Mean_Flagellum(Xtail_0, Ytail_0, Ztail_0, dir_align);   %Xtail_align: Ns*1.  
    [xup1,xup2,xup3, xdown1,xdown2,xdown3] = Generate_Virtual_Head_UpDownSides;  %xup1: M1*M2.
    Neck = [min((xup1(:)))  mean((xup2(:)))  min((xup3(:)))];
    Xtail = Xtail + Neck(1);  
    Ytail = Ytail + Neck(2);
    Ztail = Ztail + Neck(3);
    
         
    
    %================================================================================================================================================
    % Plot.
    fig = 1;
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////
    bound_height = 10; %1; %
    
    figure(fig)
    clf; set(gcf,'color','w');
    hold on; axis equal; %axis off; 
    %view([125 10]); % for full flow.
    %view([260 15]);  %view([125 10]); % for cross-sectional flow at x-planes.
    %view([-30 15]); %view([125 10]); % for cross-sectional flow at y-planes.
    view([125 10]); % view([155 20]); %for cross-sectional flow at z-planes.
    
    XL=[min(Xg(:)) max(Xg(:))];  
    %YL=[-0.4, 0.4]; ZL=[-0.4, 0.2];  % for selected full flow: no boundary.
    %YL=[-0.8, 1]; ZL=[-bound_height, 0.5];  % for selected full flow: with boundary height 0.2L, 0.5L, L.
    YL=[min(Yg(:))  max(Yg(:))]; ZL=[min(Zg(:)) max(Zg(:))];   %for cross-section flows.
    %XL=[-0.1 0.5]; YL=[-0.1 0.5]; ZL=[-0.1 0.5];   %for cross-section flow: XYZ.
    xlim([XL(1) XL(2)]); ylim([YL(1) YL(2)]); zlim([ZL(1) ZL(2)]); 
    xlabel('x'); ylabel('y'); zlabel('z');


if 0   
    % Plot sperm.
    plot3(Xtail,Ytail,Ztail,'color',[36 100 171]./255,'linewidth',2);    %'color',[0.8500    0.3250    0.0980]
    CO_up(:,:,1) = 0.6*ones(size(xup1)); CO_up(:,:,2) = 0*ones(size(xup1)); CO_up(:,:,3) = 0*ones(size(xup1));
    CO_down(:,:,1) = 0*ones(size(xdown1)); CO_down(:,:,2) = 0.6*ones(size(xdown1)); CO_down(:,:,3) = 0*ones(size(xdown1));
    s1= surf(xup1,xup2,xup3,CO_up);
    s2= surf(xdown1,xdown2,xdown3,CO_down);
    s1.EdgeColor='none' ;
    s2.EdgeColor='none' ;
    
    
    % Plot flow streamlines.
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    % Selected full flow:
    %[sx_3D,sy_3D,sz_3D] = meshgrid(xg(8:17),yg(10:12),zg(10:12)); % without boundary.            
    %[sx_3D,sy_3D,sz_3D] = meshgrid(xg(8:17),yg(10:13),zg(10:13)); % with boundary height 0.2L, 0.5L.            
    %[sx_3D,sy_3D,sz_3D] = meshgrid(xg(8:17),yg(9:14),zg(10:13)); % with boundary height L.  
    % Cross-sectional (cs) flow:
    if 1 % Full cs flows (x,y,z).
        [sx_3D,sy_3D,sz_3D] = meshgrid(xg(15),yg(5:1:end-4),zg(5:1:end-4));
        Generate_FlowStreamlines(Xg,Yg,Zg,Ug_mean,Vg_mean,Wg_mean,sx_3D,sy_3D,sz_3D, vel_dim);
        [sx_3D,sy_3D,sz_3D] = meshgrid(xg(10:1:16),yg(11),zg(5:1:end-4));
        Generate_FlowStreamlines(Xg,Yg,Zg,Ug_mean,Vg_mean,Wg_mean,sx_3D,sy_3D,sz_3D, vel_dim);
        [sx_3D,sy_3D,sz_3D] = meshgrid(xg(10:1:16),yg(5:1:end-4),zg(11));
        Generate_FlowStreamlines(Xg,Yg,Zg,Ug_mean,Vg_mean,Wg_mean,sx_3D,sy_3D,sz_3D, vel_dim);
    else % Single cs flows.
        % cs x. Two views: view([125 10]); view([260 15]);
        %[sx_3D,sy_3D,sz_3D] = meshgrid(xg(16),yg(1:1:end),zg(1:1:end)); % without boundary, vortex cross sections: xg=16/ x=0 (proximal).
        %[sx_3D,sy_3D,sz_3D] = meshgrid(xg(11),yg(1:1:end),zg(1:1:end)); % without boundary, vortex cross sections: xg=11/ x=-0.5 (middle).
        %[sx_3D,sy_3D,sz_3D] = meshgrid(xg(6),yg(1:1:end),zg(1:1:end)); % without boundary, vortex cross sections: xg=6/ x=-1 (distal).
        % cs y.
        %[sx_3D,sy_3D,sz_3D] = meshgrid(xg(1:1:end),yg(16),zg(1:1:end)); % without boundary, vortex cross sections: yg=16/ y=0.5 (left).
        %[sx_3D,sy_3D,sz_3D] = meshgrid(xg(1:1:end),yg(11),zg(1:1:end)); % without boundary, vortex cross sections: yg=11/ y=0 (middle).
        %[sx_3D,sy_3D,sz_3D] = meshgrid(xg(1:1:end),yg(6),zg(1:1:end)); % without boundary, vortex cross sections: yg=6/ y=-0.5 (right).
        % cs z.
        %[sx_3D,sy_3D,sz_3D] = meshgrid(xg(1:1:end),yg(1:1:end),zg(16)); % without boundary, vortex cross sections: zg=16/ z=0.5 (top).
        %[sx_3D,sy_3D,sz_3D] = meshgrid(xg(1:1:end),yg(1:1:end),zg(11)); % without boundary, vortex cross sections: zg=11/ z=0 (middle).
        [sx_3D,sy_3D,sz_3D] = meshgrid(xg(1:1:end),yg(1:1:end),zg(6)); % without boundary, vortex cross sections: zg=6/ z=-0.5 (bottom).
        Generate_FlowStreamlines(Xg,Yg,Zg,Ug_mean,Vg_mean,Wg_mean,sx_3D,sy_3D,sz_3D, vel_dim);
    end  
    
    % Plot the boundary surface.
    XL=[xl_field(1) xl_field(2)];
    YL=[yl_field(1) yl_field(2)];
    ZL=[-bound_height zl_field(2)]; 
    [xgrid,ygrid,zgrid] = meshgrid(linspace(XL(1),XL(2),10),...
          linspace(YL(1),YL(2),10),...
          linspace(ZL(1),ZL(2),10));                     
    CO(:,:,1) = ones(10)*.2;  %CO(:,:,1) = ones(10)*.0; % red
    CO(:,:,2) = ones(10)*.2; %CO(:,:,2) = ones(10)*.4470; % green
    CO(:,:,3) = ones(10)*0.2;  %CO(:,:,3) = ones(10)*0.7410; % blue
    surfZ = surf(xgrid(:,:,1),ygrid(:,:,1),zgrid(:,:,1), CO,'FaceAlpha',0.1);  
    surfZ.EdgeColor='none'; 
    
    
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    % Plot 2D slice planes (for no boundary case). 
    if 0 
    slice_X = 0:-0.5:-1;
    for i_slice = 1:length(slice_X)
        XL=[slice_X(i_slice) xl_field(2)];
        YL=[yl_field(1) yl_field(2)];
        ZL=[zl_field(1) zl_field(2)]; 
        [xgrid,ygrid,zgrid] = meshgrid(linspace(XL(1),XL(2),10),...
                    linspace(YL(1),YL(2),10),...
                    linspace(ZL(1),ZL(2),10));     
        CO(:,:,1) = ones(10)*.0; % red 
        CO(:,:,2) = ones(10)*.4470; % green 
        CO(:,:,3) = ones(10)*0.7410; % blue 
        surfX = surf(reshape(xgrid(:,1,:),[10 10]),reshape(ygrid(:,1,:),[10 10]),reshape(zgrid(:,1,:),[10,10]), CO,'FaceAlpha',0.1);  
        surfX.EdgeColor='none';
    end
    elseif 0
    slice_Y = -0.5:0.5:0.5;
    for i_slice = 1:length(slice_Y)
        XL=[xl_field(1) xl_field(2)];
        YL=[slice_Y(i_slice) yl_field(2)];
        ZL=[zl_field(1) zl_field(2)]; 
        [xgrid,ygrid,zgrid] = meshgrid(linspace(XL(1),XL(2),10),...
                    linspace(YL(1),YL(2),10),...
                    linspace(ZL(1),ZL(2),10));     
        CO(:,:,1) = ones(10)*.0; % red 
        CO(:,:,2) = ones(10)*.4470; % green 
        CO(:,:,3) = ones(10)*0.7410; % blue 
        surfY = surf(reshape(xgrid(1,:,:),[10 10]),reshape(ygrid(1,:,:),[10 10]),reshape(zgrid(1,:,:),[10,10]), CO,'FaceAlpha',0.1);  
        surfY.EdgeColor='none';
    end
    elseif 1
    slice_Z = -0.5:0.5:0.5;
    for i_slice = 1:length(slice_Z)
        XL=[xl_field(1) xl_field(2)];
        YL=[yl_field(1) yl_field(2)];
        ZL=[slice_Z(i_slice) zl_field(2)]; 
        [xgrid,ygrid,zgrid] = meshgrid(linspace(XL(1),XL(2),10),...
                    linspace(YL(1),YL(2),10),...
                    linspace(ZL(1),ZL(2),10));     
        CO(:,:,1) = ones(10)*.0; % red 
        CO(:,:,2) = ones(10)*.4470; % green 
        CO(:,:,3) = ones(10)*0.7410; % blue 
        surfZ = surf(reshape(xgrid(:,:,1),[10 10]),reshape(ygrid(:,:,1),[10 10]),reshape(zgrid(:,:,1),[10,10]), CO,'FaceAlpha',0.1);  
        surfZ.EdgeColor='none';
    end
    end
    
    shading interp
    colorbar; 
    %camlight right
    %lighting gouraud
     
    
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    caxis([0 1]); % caxis([0 5]); 
    %return
    
end    
    if 1
    % Plot Cartesian axes.
    %O = [0.2 0.2 -0.3]; %for full flow: no boundary case.////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    %O = [0.35 -0.6 -0.15]; %for full flow: boundary height 0.2L.
    %O = [-0.5 0 -0.2]; %for full flow: boundary height 0.5L, L.
    O = [-1.3 -0.8 -0.8]; %for cross section: proximal, middle, distal.
    Cx=[1; 0; 0]; 
    Cy=[0; 1; 0]; 
    Cz=[0; 0; 1]; 
    C=[Cx Cy Cz];
    Plot_Cartesian_Axes(fig,O,C)  
    end
    
    camlight right
    lighting gouraud
     
    
    return
end




%% Fig. 2: streamslice.

if 0
    
    %=====================================================================================================
    % Get flow field.
    load NearFlow_CF_sp8_NoBound.mat;
    
    xg=linspace(xl_field(1),xl_field(2),Nx_field);
    yg=linspace(yl_field(1),yl_field(2),Ny_field); 
    zg=linspace(zl_field(1),zl_field(2),Nz_field); 
    [Xg,Yg,Zg]=meshgrid(xg,yg,zg);   
    
    Ug_mean = mean(Ug,2);
    Vg_mean = mean(Vg,2);
    Wg_mean = mean(Wg,2);
    Ug_mean = reshape(Ug_mean,size(Xg,1),size(Xg,2),size(Xg,3));
    Vg_mean = reshape(Vg_mean,size(Xg,1),size(Xg,2),size(Xg,3));
    Wg_mean = reshape(Wg_mean,size(Xg,1),size(Xg,2),size(Xg,3));
    
    vel_raw=(Ug_mean.*Ug_mean+Vg_mean.*Vg_mean+Wg_mean.*Wg_mean).^0.5; 
    vel_nd=vel_raw.*2*pi;  %in units of 'flagellum length/ beat cycle',
    vel_dim = vel_nd.*scale_vel;  % in units of '\mu m/s'.
   
    
    
    %=====================================================================================================
    % Get sperm shape.
    load XNodes_sp8_NoBoundStokeslet.mat;
    
    nt_XNodes=find(t_nd==t_nd_flow(end));
    Q=Nhh+Ns;   
    Xtail_0 = Xs(Nhh+1:Q,1:nt_XNodes); 
    Ytail_0 = Xs(Q+Nhh+1:2*Q,1:nt_XNodes); 
    Ztail_0 = Xs(2*Q+Nhh+1:3*Q,1:nt_XNodes); 
    
    % Get average sperm body shape.
    dir_align = [1, 0, 0];
    [Xtail, Ytail, Ztail] = Get_Aligned_Mean_Flagellum(Xtail_0, Ytail_0, Ztail_0, dir_align);   %Xtail_align: Ns*1.  
    [xup1,xup2,xup3, xdown1,xdown2,xdown3] = Generate_Virtual_Head_UpDownSides;  %xup1: M1*M2.
    Neck = [min((xup1(:)))  mean((xup2(:)))  min((xup3(:)))];
    Xtail = Xtail + Neck(1);  
    Ytail = Ytail + Neck(2);
    Ztail = Ztail + Neck(3);
    
         
    
    %================================================================================================================================================
    % Plot.
    fig=3;
    bound_height = 10; %0.2; 
    
    figure(fig)
    clf; set(gcf,'color','w');
    hold on; axis equal; axis off; 
    view([90 0]); 
    XL=[min(Xg(:)) max(Xg(:))]; YL=[min(Yg(:))-0.07  max(Yg(:))]; ZL=[max([-bound_height min(Zg(:))])-0.07  max(Zg(:))];  
    xlim([XL(1) XL(2)]); ylim([YL(1) YL(2)]); zlim([ZL(1) ZL(2)]); 
    
    if 1
    % Plot sperm.
    plot3(Xtail,Ytail,Ztail,'color',[36 100 171]./255,'linewidth',2);    %'color',[0.8500    0.3250    0.0980]
    CO_up(:,:,1) = 0.6*ones(size(xup1)); CO_up(:,:,2) = 0*ones(size(xup1)); CO_up(:,:,3) = 0*ones(size(xup1));
    CO_down(:,:,1) = 0*ones(size(xdown1)); CO_down(:,:,2) = 0.6*ones(size(xdown1)); CO_down(:,:,3) = 0*ones(size(xdown1));
    s1= surf(xup1,xup2,xup3,CO_up);
    s2= surf(xdown1,xdown2,xdown3,CO_down);
    s1.EdgeColor='none' ;
    s2.EdgeColor='none' ;
    
    
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    % Plot cross-sectional streamlines (x=0,-0.5,-1).
    slice(Xg,Yg,Zg,vel_dim,[-1],[],[]);
    ss=streamslice(Xg,Yg,Zg,6*Ug_mean,6*Vg_mean,6*Wg_mean,[-1],[],[],0.2); 
    set(ss,'color','w','linewidth',15);
    
    shading interp;
    colormap(gca,'parula')
    colorbar; 
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    caxis([0 0.3])
    return
    end
    % Plot Cartesian axes.
    O = [0 -1 -1]; %for cross section.
    Cx=[1; 0; 0]; 
    Cy=[0; 1; 0]; 
    Cz=[0; 0; 1]; 
    C=[Cx Cy Cz];
    Plot_Cartesian_Axes(fig,O,C)  
    
    %camlight right
    %lighting gouraud
  
    
    return
end




    



%% Fig. 3: Colorbar.

figure(3)
clf; set(gcf,'color','w');
hold on; axis equal; axis off; box off;
    
colormap(parula); 
caxis([0,1]);

c = colorbar('Ticks',[0,1],'TickLabels',{'0','1'},'FontSize',50,'fontname','Times');  %,'TickLabels',{'0','1'}
c.Label.String = 'Flow velocity (\mum/s)';
c.TickLabelInterpreter='latex';
c.Position = [0.7 0.1 0.01 0.65];  % orientation: vertical
%c.Position = [0.85 0.2 0.01 0.7];  % orientation: vertical

return

c = colorbar('southoutside','Ticks',[0,1],'TickLabels',{'0','1'},'FontSize',80,'FontName','Times');  
c.Label.String = 'Flow velocity';
c.TickLabelInterpreter='latex'; 
c.Position = [0.2 0.3 0.45 0.015];  % orientation: horizontal


c = colorbar('southoutside','Ticks',[0,1],'TickLabels',{'0','1'},'FontSize',80,'FontName','Times');  
c.Label.String = '(\mum/s)';
c.TickLabelInterpreter='latex'; 
c.Position = [0.2 0.5 0.45 0.015];  % orientation: horizontal

