function rxm_IdealFlow_AverageNear_Plot

% This function plots the average flow field around the reconstructed human sperm, relative to the comoving frame of reference.
% Figures:
% 1. 3D average flow. 
% 2. 2D slicestream.
% 3. Colorbar.

clc




%% Fig. 1: 3D flow.

if 1
    %=====================================================================================================
    % Get flow field.
    load NearFlow_CF_alpha0.6.mat;
    
    xg=linspace(xl_field(1),xl_field(2),Nx_field);
    yg=linspace(yl_field(1),yl_field(2),Ny_field); 
    zg=linspace(zl_field(1),zl_field(2),Nz_field); 
    [Xg,Yg,Zg]=meshgrid(xg,yg,zg);   
    
    Ug_mean = mean(Ug,2);  %Ug_mean = mean(Ug(:,1:20),2); %
    Vg_mean = mean(Vg,2);  %Vg_mean = mean(Vg(:,1:20),2);  %
    Wg_mean = mean(Wg,2);  %Wg_mean = mean(Wg(:,1:20),2);  %
    Ug_mean = reshape(Ug_mean,size(Xg,1),size(Xg,2),size(Xg,3));
    Vg_mean = reshape(Vg_mean,size(Xg,1),size(Xg,2),size(Xg,3));
    Wg_mean = reshape(Wg_mean,size(Xg,1),size(Xg,2),size(Xg,3));
    
    vel_raw=(Ug_mean.*Ug_mean+Vg_mean.*Vg_mean+Wg_mean.*Wg_mean).^0.5; 
    vel_nd=vel_raw.*2*pi;  %in units of 'flagellum length/ beat cycle',
    
    
    
    %=====================================================================================================
    % Get sperm shape.
    load XNodes_Lab_alpha0.6.mat;
    nt_XNodes=find(tRange==tRange_flow(end));
    Q=Nhh+Ns;   
    Xtail_0 = Xs(Nhh+1:Q,1:nt_XNodes); 
    Ytail_0 = Xs(Q+Nhh+1:2*Q,1:nt_XNodes); 
    Ztail_0 = Xs(2*Q+Nhh+1:3*Q,1:nt_XNodes); 
    
    % Get average sperm body shape.
    dir_align = [-1, 0, 0];
    [Xtail, Ytail, Ztail] = Get_Aligned_Mean_Flagellum(Xtail_0, Ytail_0, Ztail_0, dir_align);   %Xtail_align: Ns*1.  
    [xup1,xup2,xup3, xdown1,xdown2,xdown3] = Generate_Virtual_Head_UpDownSides;  %xup1: M1*M2.
    Neck = [min((xup1(:)))  mean((xup2(:)))  min((xup3(:)))];
    Xtail = Xtail + Neck(1);  
    Ytail = Ytail + Neck(2);
    Ztail = Ztail + Neck(3);
    
         
    
    %================================================================================================================================================
    % Plot.
    fig = 2;
    
    figure(fig)
    clf; set(gcf,'color','w');
    hold on; axis equal; axis off; 
    xlabel('x'); ylabel('y'); zlabel('z');
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////
    XL=[xl_field(1) xl_field(2)]; YL=[-0.6 0.6]; ZL=[-0.5 0.5]; % for full flow (alpha 0.6). 
    %  XL=[-0.2 xl_field(2)]; YL=[min(Yg(:))  max(Yg(:))]; ZL=[min(Zg(:)) max(Zg(:))];   %for cross-section flow (alpha 0.6): xg=15 (x=0.9, distal).
    % XL=[xl_field(1) xl_field(2)]; YL=[min(Yg(:))  max(Yg(:))]; ZL=[min(Zg(:)) max(Zg(:))];   %for cross-section flow (alpha 0.6): xg=11 (x=0.5, middle).
    % XL=[xl_field(1) xl_field(2)]; YL=[min(Yg(:))  max(Yg(:))]; ZL=[min(Zg(:)) max(Zg(:))];   %for cross-section flow (alpha 0.6): xg=6 (x=0, proximal).    
    xlim([XL(1) XL(2)]); ylim([YL(1) YL(2)]); zlim([ZL(1) ZL(2)]); 
    view([-45 20]);

    
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
    if 1 % Plot 2D slice planes (for full flow).
    slice_X = [0 0.5 0.9]; 
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
    end
    
    
    % Plot flow streamlines.
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////
    [sx_3D,sy_3D,sz_3D] = meshgrid(xg(7:10),yg(7:15),zg(10:12));  % for full flow (alpha 0.6).
    %[sx_3D,sy_3D,sz_3D] = meshgrid(xg(6),yg(1:1:end),zg(1:1:end)); % for vortex cross sections: xg=6 (x=0, proximal).
    %[sx_3D,sy_3D,sz_3D] = meshgrid(xg(11),yg(1:1:end),zg(1:1:end)); % for vortex cross sections: xg=11 (x=0.5, middle).
    %[sx_3D,sy_3D,sz_3D] = meshgrid(xg(15),yg(1:1:end),zg(1:1:end)); % for vortex cross sections: xg=15 (x=0.9, distal).
    Generate_FlowStreamlines(Xg,Yg,Zg,Ug_mean,Vg_mean,Wg_mean,sx_3D,sy_3D,sz_3D, vel_nd);
    shading interp
    
    
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////
    % Arrows indicating the flow directions.
    if 1   
        [sx_cone1,sy_cone1,sz_cone1] = meshgrid(xg(4),yg(8:3:14),zg([8 11 15]));  
        sx_cone1([1 3],:,[1 3])=NaN; sy_cone1([1 3],:,[1 3])=NaN; sz_cone1([1 3],:,[1 3])=NaN;
        sx_cone1(2,:,2)=NaN; sy_cone1(2,:,2)=NaN; sz_cone1(2,:,2)=NaN;
        [sx_cone2,sy_cone2,sz_cone2] = meshgrid(xg(13),yg(8:3:14),zg([8 11 15])); 
        sx_cone2([1 3],:,[1 3])=NaN; sy_cone2([1 3],:,[1 3])=NaN; sz_cone2([1 3],:,[1 3])=NaN;
        sx_cone2(2,:,2)=NaN; sy_cone2(2,:,2)=NaN; sz_cone2(2,:,2)=NaN;
        scal1 = 15; scal2 = 15;   
        Plot_StreamCone(fig,Xg,Yg,Zg,Ug_mean,Vg_mean,Wg_mean,sx_cone1,sy_cone1,sz_cone1,scal1)
        Plot_StreamCone(fig,Xg,Yg,Zg,Ug_mean,Vg_mean,Wg_mean,sx_cone2,sy_cone2,sz_cone2,scal2)
    end
    caxis([0 0.005]);  % for full flow (alpha 0.6).
    %caxis([0 0.001]); % for vortex cross sections.
    colorbar; 
    
    camlight right
    lighting gouraud
    return
end

    % Plot Cartesian axes.
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////
    O = [0 -0.4 -0.4]; % for full flow (alpha 0.6).
    %O = [-0.6 0.5 -0.95]; % for vortex cross sections: xg=6 (x=0, proximal).
    %O = [0 -0.2 -0.95]; %for cross section: xg=11 (x=0.5, middle).
    %O = [0.3 0.6 -0.8]; %for cross section: xg=15 (x=0.9, distal).
    Cx=[-1; 0; 0]; 
    Cy=[0; -1; 0]; 
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
    load NearFlow_CF_alpha0.6.mat;
    
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
    
    
    %=====================================================================================================
    % Get sperm shape.
    load XNodes_Lab_alpha0.6.mat;
    nt_XNodes=find(tRange==tRange_flow(end));
    Q=Nhh+Ns;   
    Xtail_0 = Xs(Nhh+1:Q,1:nt_XNodes); 
    Ytail_0 = Xs(Q+Nhh+1:2*Q,1:nt_XNodes); 
    Ztail_0 = Xs(2*Q+Nhh+1:3*Q,1:nt_XNodes); 
    
    % Get average sperm body shape.
    dir_align = [-1, 0, 0];
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
    xlabel('x'); ylabel('y'); zlabel('z');
    view([-90 0]); 
    XL=[min(Xg(:)) max(Xg(:))]; YL=[min(Yg(:))  max(Yg(:))]; ZL=[min(Zg(:))  max(Zg(:))]; 
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
    
    
    % Plot cross-sectional streamlines (x=0,0.5,0.9).
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    slice(Xg,Yg,Zg,vel,[0.9],[],[]);
    %ss=streamslice(Xg,Yg,Zg,6*Ug_mean,6*Vg_mean,6*Wg_mean,[0.9],[],[],0.3); 
    %set(ss,'color','w','linewidth',10);
    ss=streamslice(Xg,Yg,Zg,6*Ug_mean,6*Vg_mean,6*Wg_mean,[0.9],[],[],0.2); 
    set(ss,'color','w','linewidth',15);
    shading interp;
    colormap(gca,'parula')
    colorbar; 
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    %caxis([0 0.005]);  % for full flow (alpha 0.6).
    %caxis([0 0.001]); % for vortex cross sections: xg=6 (x=0, proximal) and xg=11 (x=0.5, middle).
    caxis([0 0.0005]); % for cross-sectional flow: xg=15, x=0.9 (distal).
else
    % Plot Cartesian axes.
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    O = [0 1 -1]; %for cross section: proximal.
    Cx=[-1; 0; 0]; 
    Cy=[0; -1; 0]; 
    Cz=[0; 0; 1]; 
    C=[Cx Cy Cz];
    Plot_Cartesian_Axes(fig,O,C)      
end

return
end




    



%% Fig. 3: Colorbar.

figure(3)
clf; set(gcf,'color','w');
hold on; axis equal; axis off; box off;
       
colormap(parula); 
clim = 0.005;
caxis([0,clim]);

c = colorbar('Ticks',[0,clim],'FontSize',60,'TickLabels',{'0','0.005'},'fontname','Times');  %,'TickLabels',{'0','1'}
c.Label.String = 'Flow velocity';
c.TickLabelInterpreter='latex';
c.Position = [0.7 0.2 0.01 0.7];  % orientation: vertical


c = colorbar('southoutside','Ticks',[0,clim],'TickLabels',{'0','0.001'},'FontSize',60,'FontName','Times');  
c.Label.String = 'Flow velocity';
c.TickLabelInterpreter='latex'; 
c.Position = [0.1 0.3 0.55 0.015];  % orientation: horizontal


