function rxm_ReconFlow_AverageFar_Plot

% This function plots the far-field average flow around the reconstructed human sperm, relative to the comoving frame of reference.
% Figures:
% 1. 3D average flow. 
% 2. 2D slicestream.
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
    load FarFlow_CF_sp8_NoBound.mat;
    
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
    
    figure(fig)
    clf; set(gcf,'color','w');
    hold on; axis equal; axis off; 
    %//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    XL=[xl_field(1) xl_field(2)]; YL=[yl_field(1) yl_field(2)]; ZL=[zl_field(1) zl_field(2)];  % for flow.
    %XL=[0.7 1.5]; YL=[-0.2 0.5]; ZL=[-0.2 0.5];   % for XYZ axes. 
    xlim([XL(1) XL(2)]); ylim([YL(1) YL(2)]); zlim([ZL(1) ZL(2)]); 
    view([155 15]);  
    xlabel('x'); ylabel('y'); zlabel('z');
    
if 1    
    % Plot sperm.
    plot3(Xtail,Ytail,Ztail,'color',[36 100 171]./255,'linewidth',2);    %'color',[0.8500    0.3250    0.0980]
    CO_up(:,:,1) = 0.6*ones(size(xup1)); CO_up(:,:,2) = 0*ones(size(xup1)); CO_up(:,:,3) = 0*ones(size(xup1));
    CO_down(:,:,1) = 0*ones(size(xdown1)); CO_down(:,:,2) = 0.6*ones(size(xdown1)); CO_down(:,:,3) = 0*ones(size(xdown1));
    s1= surf(xup1,xup2,xup3,CO_up);
    s2= surf(xdown1,xdown2,xdown3,CO_down);
    s1.EdgeColor='none' ;
    s2.EdgeColor='none' ;
    
    
    if 1 % Plot 2D slice planes.
    % X planes.
    XL1_mat = [-4 -1.5 0 1 4];
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
    end
    
    
    % Plot full flow streamlines.
        % cross-sectional flow at X=0. 
        theta = linspace(0,2*pi,61);
        sy = 1*cos(theta);  
        sz = 1*sin(theta);
        sx = zeros(size(sy))+0;
        Generate_FlowStreamlines(Xg,Yg,Zg,Ug_mean,Vg_mean,Wg_mean,sx,sy,sz, vel_dim);       
        % cross-sectional flow at Y=0.
        r = linspace(-1,1,41); 
        nr = length(r);
        xr = -10:0.5:9; 
        nx = length(xr);
        sx = nan(nr,nx); 
        for i_x = 1:nx
           sx(:,i_x) = xr(i_x)*r';
        end
        sz = repmat( linspace(-10,10,nr)',1,nx ); 
        sy = zeros(size(sx))+0;
        Generate_FlowStreamlines(Xg,Yg,Zg,Ug_mean,Vg_mean,Wg_mean,sx,sy,sz, vel_dim);
        % cross-sectional flow at Z=0.
        r = linspace(-1,1,41); 
        nr = length(r);
        xr = -8:0.5:10;   
        nx = length(xr);
        sx = nan(nr,nx); 
        for i_x = 1:nx
           sx(:,i_x) = xr(i_x)*r';
        end
        sy = repmat( linspace(-10,10,nr)',1,nx ); 
        sz = zeros(size(sx))+0;
        Generate_FlowStreamlines(Xg,Yg,Zg,Ug_mean,Vg_mean,Wg_mean,sx,sy,sz, vel_dim); 
         
        shading interp
        colorbar; 
        caxis([0 0.01]);       
        
        
    
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////
    % Arrows indicating the streamlines' directions.
    if 1   % for full flow (reconstructed sperm, sp#8, no boundary).
        % arrows located at the top and bottom.
        %[sx_cone1,sy_cone1,sz_cone1] = meshgrid(xg([23]),yg([3 21 40]),zg([3 18 39]));  
        [sx_cone1,sy_cone1,sz_cone1] = meshgrid(xg([23]),yg([3 22 40]),zg([3 18 39]));  
        sx_cone1([1 3],:,:)=NaN; sy_cone1([1 3],:,:)=NaN; sz_cone1([1 3],:,:)=NaN;
        sx_cone1(2,:,2)=NaN; sy_cone1(2,:,2)=NaN; sz_cone1(2,:,2)=NaN;
        % arrows located at the left.
        [sx_cone2,sy_cone2,sz_cone2] = meshgrid(xg([10 28 40]),yg([2 26 40]),zg([22]));  
        sx_cone2(:,[1 3],:)=NaN; sy_cone2(:,[1 3],:)=NaN; sz_cone2(:,[1 3],:)=NaN;
        sx_cone2([1 2],2,:)=NaN; sy_cone2([1 2],2,:)=NaN; sz_cone2([1 2],2,:)=NaN;
        % arrows located in the front.
        %[sx_cone3,sy_cone3,sz_cone3] = meshgrid(xg([6 23 36]),yg([22]),zg([8 22 40]));  
        [sx_cone3,sy_cone3,sz_cone3] = meshgrid(xg([6 23 36]),yg([22]),zg([8 24 40]));  
        sx_cone3(:,:,[1 3])=NaN;  sy_cone3(:,:,[1 3])=NaN;   sz_cone3(:,:,[1 3])=NaN;
        sx_cone3(:,[1 2],2)=NaN;  sy_cone3(:,[1 2],2)=NaN;   sz_cone3(:,[1 2],2)=NaN;
        % arrows located at the back.
        %[sx_cone4,sy_cone4,sz_cone4] = meshgrid(xg([6 23 36]),yg([25]),zg([8 24 40]));  
        [sx_cone4,sy_cone4,sz_cone4] = meshgrid(xg([6 23 36]),yg([25]),zg([8 23 40]));  
        sx_cone4(:,:,[1 3])=NaN;  sy_cone4(:,:,[1 3])=NaN;   sz_cone4(:,:,[1 3])=NaN;
        sx_cone4(:,[2 3],2)=NaN;  sy_cone4(:,[2 3],2)=NaN;   sz_cone4(:,[2 3],2)=NaN;
        
        scal1 = 15000; scal2 = 25000; scal3 = 4500; scal4 = 8000;   
        Plot_StreamCone(fig,Xg,Yg,Zg,Ug_mean,Vg_mean,Wg_mean,sx_cone1,sy_cone1,sz_cone1,scal1)
        Plot_StreamCone(fig,Xg,Yg,Zg,Ug_mean,Vg_mean,Wg_mean,sx_cone2,sy_cone2,sz_cone2,scal2)
        Plot_StreamCone(fig,Xg,Yg,Zg,Ug_mean,Vg_mean,Wg_mean,sx_cone3,sy_cone3,sz_cone3,scal3)
        Plot_StreamCone(fig,Xg,Yg,Zg,Ug_mean,Vg_mean,Wg_mean,sx_cone4,sy_cone4,sz_cone4,scal4)
    end
    
    camlight right
    lighting gouraud
    return
end


% Plot Cartesian axes.
    % Plot Cartesian axes.
    O =[1 0 0]; % for PCA modes.
    Cx=[1; 0; 0]; 
    Cy=[0; 1; 0]; 
    Cz=[0; 0; 1]; 
    C=[Cx Cy Cz];
    Plot_Cartesian_Axes(fig,O,C)  
    camlight right
    lighting gouraud    
    return
end




%% Fig. 2: streamslice.

if 0
    
    %=====================================================================================================
    % Get flow field.
    load FarFlow_CF_sp8_NoBound.mat;   
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
    vel=vel_raw.*2*pi;  %in units of 'flagellum length/ beat cycle',
    vel_dim = vel.*scale_vel;  % in units of '\mu m/s'.
   
    
    
    %=====================================================================================================
    % Get sperm shape.
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
    fig = 4;
    
    figure(fig)
    clf; set(gcf,'color','w');
    hold on; axis equal; axis off; 
    %view([155 15]);   
    view([90 0]); % for YZ-plane view
    %view([0 0]); % for XZ-plane view
    %view([0 90]); % for XY-plane view
    XL_temp=[min(Xg(:)) max(Xg(:))]; YL_temp=[min(Yg(:))  max(Yg(:))]; ZL_temp=[min(Zg(:))  max(Zg(:))]; 
    xlim([XL_temp(1) XL_temp(2)]); ylim([YL_temp(1) YL_temp(2)]); zlim([ZL_temp(1) ZL_temp(2)]); 
   
    % Plot sperm.
    plot3(Xtail,Ytail,Ztail,'color',[36 100 171]./255,'linewidth',2);    %'color',[0.8500    0.3250    0.0980]
    CO_up(:,:,1) = 0.6*ones(size(xup1)); CO_up(:,:,2) = 0*ones(size(xup1)); CO_up(:,:,3) = 0*ones(size(xup1));
    CO_down(:,:,1) = 0*ones(size(xdown1)); CO_down(:,:,2) = 0.6*ones(size(xdown1)); CO_down(:,:,3) = 0*ones(size(xdown1));
    s1= surf(xup1,xup2,xup3,CO_up);
    s2= surf(xdown1,xdown2,xdown3,CO_down);
    s1.EdgeColor='none' ;
    s2.EdgeColor='none' ;
    
    
    % Plot cross-sectional streamlines (x=-4,-1.5,0,1,4; y=0; z=0).
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    slice(Xg,Yg,Zg,vel_dim,[4],[],[]);
    ss=streamslice(Xg,Yg,Zg,6*Ug_mean,6*Vg_mean,6*Wg_mean,[4],[],[],0.3); %0.2 %0.1 %0.3
    set(ss,'color','w','linewidth',10); %20 %10
    shading interp;
    colormap(gca,'parula')
    colorbar; 
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    caxis([0 0.1])
    return
    
    % Plot Cartesian axes.
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    O = [0 1 -1]; %for cross section: proximal.
    Cx=[-1; 0; 0]; 
    Cy=[0; -1; 0]; 
    Cz=[0; 0; 1]; 
    C=[Cx Cy Cz];
    Plot_Cartesian_Axes(fig,O,C)  
    
    
    return
end




    



%% Fig. 3: Colorbar.

figure(3)
clf; set(gcf,'color','w');
hold on; axis equal; axis off; box off;
       
colormap(parula); 
clim = 0.005;
caxis([0,clim]);

c = colorbar('Ticks',[0,clim],'FontSize',30,'TickLabels',{'0','0.1'},'fontname','Times'); %,'FontSize',65 %'TickLabels',{'0','$5\times10^{-5}$'}
c.Label.String = 'Flow velocity (\mum/s)';
c.TickLabelInterpreter='latex';
c.Position = [0.7 0.1 0.005 0.8];  % orientation: vertical
%c.Position = [0.7 0.1 0.01 0.5];  % orientation: vertical

return
c = colorbar('southoutside','Ticks',[0,clim],'TickLabels',{'0','$5\times10^{-5}$'},'FontSize',80,'FontName','Times');  
c.Label.String = 'Flow velocity';
c.TickLabelInterpreter='latex'; 
c.Position = [0.1 0.3 0.44 0.015];  % orientation: horizontal


