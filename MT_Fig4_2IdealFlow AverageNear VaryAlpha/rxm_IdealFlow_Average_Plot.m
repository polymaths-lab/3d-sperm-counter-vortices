function rxm_IdealFlow_Average_Plot

% This function plots the average flow field around the reconstructed human sperm, relative to the comoving frame of reference.
% Figures:
% 1. 3D average flow. 
% 2. Colorbar.

clc




%% Fig. 1: 3D flow.

if 1
    %=====================================================================================================
    % Get flow field.
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////
    load NearFlow_CF_alpha0.4.mat;
    
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
    
    
    
    %=====================================================================================================
    % Get sperm shape.
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////
    load XNodes_Lab_alpha0.4.mat;
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
    fig=1;
    figure(fig)
    clf; set(gcf,'color','w');
    hold on; axis equal; axis off;   
    XL=[xl_field(1) xl_field(2)]; 
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////
    %view([-45 20]); YL=[-0.6 0.6]; ZL=[-0.5 0.5]; % for alpha 0.6. 
    %view([-45 20]); YL=[-0.6 0.6]; ZL=[-0.7 0.7]; % for alpha1.     
    %view([-40 40]); YL=[-0.6 0.6]; ZL=[-0.4 0.4]; % for alpha0.     
    view([-45 20]); YL=[-0.5 0.5]; ZL=[-0.45 0.55]; % for alpha0.4.  
    %YL=[-0.6 0.6]; ZL=[-0.7 0.7]; % for XYZ. 
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
    
    % Plot flow streamlines and cones.
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////
    %[sx_3D,sy_3D,sz_3D] = meshgrid(xg(7:10),yg(7:15),zg(10:12));  % for alpha 0.6, 1.
    %[sx_3D,sy_3D,sz_3D] = meshgrid(xg(6:11),yg(7:15),zg(10:12));  % for alpha0.   
    [sx_3D,sy_3D,sz_3D] = meshgrid(xg(7:10),yg(6:16),zg(9:13));  % for alpha0.4.   
    Generate_FlowStreamlines(Xg,Yg,Zg,Ug_mean,Vg_mean,Wg_mean,sx_3D,sy_3D,sz_3D, vel_nd);
    shading interp
    colormap(gca,'parula')
    colorbar; 
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////
    %caxis([0 0.005]); [sx_cone1,sy_cone1,sz_cone1] = meshgrid(xg(4),yg(7:4:15),zg(7:4:15));  [sx_cone2,sy_cone2,sz_cone2] = meshgrid(xg(13),yg(7:4:15),zg(7:4:15)); scal1 = 20; scal2 = 20;   % for alpha 0.6.
    %caxis([0 0.005]); [sx_cone1,sy_cone1,sz_cone1] = meshgrid(xg(4),yg(7:4:15),zg(7:4:15));  [sx_cone2,sy_cone2,sz_cone2] = meshgrid(xg(13),yg(7:4:15),zg(7:4:15)); scal1 = 30; scal2 = 30;   % for alpha 1.
    %caxis([0 0.01]); [sx_cone1,sy_cone1,sz_cone1] = meshgrid(xg(4),yg(8:3:14),zg(8:3:14));  [sx_cone2,sy_cone2,sz_cone2] = meshgrid(xg(16),yg(8:3:14),zg(8:3:14));  scal1 = 20; scal2 = 20; % for alpha0.
    caxis([0 0.005]); [sx_cone1,sy_cone1,sz_cone1] = meshgrid(xg(5),yg(8:3:14),zg(8:3:14));  [sx_cone2,sy_cone2,sz_cone2] = meshgrid(xg(13),yg(8:3:14),zg(8:3:14));  scal1 = 20; scal2 = 20;  % for alpha0.4.    
    Plot_StreamCone(fig,Xg,Yg,Zg,Ug_mean,Vg_mean,Wg_mean,sx_cone1,sy_cone1,sz_cone1,scal1)
    Plot_StreamCone(fig,Xg,Yg,Zg,Ug_mean,Vg_mean,Wg_mean,sx_cone2,sy_cone2,sz_cone2,scal2)

else   
    % Plot Cartesian axes.
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////
    O = [1 0.5 -0.5]; % uniform for different alpha values (0, 0.4, 1...).
    %O = [0 -0.4 -0.4]; % for alpha 0.6, 1.
    %O = [0 -0.4 -0.35]; % for alpha0.
    %O = [0.1 -0.25 -0.4]; % for alpha0.4.
    Cx=[-1; 0; 0]; 
    Cy=[0; -1; 0]; 
    Cz=[0; 0; 1]; 
    C=[Cx Cy Cz];
    Plot_Cartesian_Axes(fig,O,C)  
    
end  
    camlight right
    lighting gouraud
     
    
    return
end




%% Fig. 2: Colorbar.

figure(2)
clf; set(gcf,'color','w');
hold on; axis equal; axis off; box off;
    
colormap(parula); 
clim = 0.005;
caxis([0,clim]);

c = colorbar('Ticks',[0,clim],'FontSize',50,'TickLabels',{'0','0.01'},'fontname','Times');  %,'TickLabels',{'0','1'}
c.Label.String = 'Flow velocity';
c.TickLabelInterpreter='latex';
c.Position = [0.85 0.2 0.01 0.7];  % orientation: vertical


c = colorbar('southoutside','Ticks',[0,clim],'TickLabels',{'0','0.01'},'FontSize',50,'FontName','Times');  
c.Label.String = 'Flow velocity';
c.TickLabelInterpreter='latex'; 
c.Position = [0.2 0.3 0.45 0.015];  % orientation: horizontal


