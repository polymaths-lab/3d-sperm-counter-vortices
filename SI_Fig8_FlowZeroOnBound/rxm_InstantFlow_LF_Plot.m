function rxm_InstantFlow_LF_Plot

% This function plots the instant flow around the reconstructed sperm, in the stationary lab frame. 
% This function aims to check our implementation correctness of the image system (¡¯Blakelet') by confirming the zero flow velocity on the boundary surface. 

% Fig 1: colorbar. 
% Fig 2: 3D streamlines + 2D slices.
% Fig 3: XYZ coordinate system.

clc



%% Fig. 1: Colorbar.
if 0

figure(1)
clf; set(gcf,'color','w');
hold on; axis equal; axis off; box off;
    
colormap(parula); 
caxis([0,1]);

c = colorbar('Ticks',[0,1],'TickLabels',{'0','20'},'FontSize',50,'fontname','Times');  %,'TickLabels',{'0','1'}
c.Label.String = 'Flow velocity (\mum/s)';
c.TickLabelInterpreter='latex';
c.Position = [0.85 0.2 0.01 0.7];  % orientation: vertical

return

c = colorbar('southoutside','Ticks',[0,1],'FontSize',80,'FontName','Times');  
c.Label.String = 'Flow velocity (\mum/s)';
c.TickLabelInterpreter='latex'; 
c.Position = [0.2 0.3 0.45 0.015];  % orientation: horizontal


return 

end


%% For figure 2.

if 0
sp=8;

load 01_lab_frame_raw_traces_2017_2018_rxmSSppxyz.mat; 
arc_mat = arclength(X{sp},Y{sp},Z{sp});
arc_mean = mean(arc_mat(end,:));

load FreeSperm_Frequency.mat;  
scale_vel = arc_mean*HF_Freq(sp); % in units o '\mu m/s'.

clearvars -except scale_vel



%% Load and process data.

% select a time point to present the instant flow.
nt0 = 15; 


load NearFlow_LF_sp8_BoundH0.2.mat; 

xg=linspace(xl_field(1),xl_field(2),Nx_field);
yg=linspace(yl_field(1),yl_field(2),Ny_field); 
zg=linspace(zl_field(1),zl_field(2),Nz_field); 
[Xg,Yg,Zg]=meshgrid(xg,yg,zg);

Ug0 = Ug(:,nt0);
Vg0 = Vg(:,nt0);
Wg0 = Wg(:,nt0);
Ug0 = reshape(Ug0,size(Xg,1),size(Xg,2),size(Xg,3));
Vg0 = reshape(Vg0,size(Xg,1),size(Xg,2),size(Xg,3));
Wg0 = reshape(Wg0,size(Xg,1),size(Xg,2),size(Xg,3));
    
vel_raw=(Ug0.*Ug0+Vg0.*Vg0+Wg0.*Wg0).^0.5; 
vel_nd=vel_raw.*2*pi;  %in units of 'flagellum length/ beat cycle',
vel_dim = vel_nd.*scale_vel;  % in units of '\mu m/s'.
   
    


load XNodes_sp8_BoundBlakelet_H0.2.mat;
nt0_X = find(t_nd==t_nd_flow(nt0));
t_nd_flow(nt0)
t_nd(nt0_X)

Q=Nhh+Ns;   
Xtail = Xs(Nhh+1:Q,nt0_X); 
Ytail = Xs(Q+Nhh+1:2*Q,nt0_X); 
Ztail = Xs(2*Q+Nhh+1:3*Q,nt0_X);
[Xtail,Ytail,Ztail] = Get_Smooth_Flag(Xtail,Ytail,Ztail);

neck_nt0 = [Xtail(1),Ytail(1),Ztail(1)];
%head_x0_nt0 = head_x0(nt0_X,:);
head_tangent_nt0=head_tangent(nt0_X,:)';
head_normal_nt0=head_normal(nt0_X,:)';
head_binormal_nt0=cross(head_tangent_nt0,head_normal_nt0);
B_nt0=[head_tangent_nt0,head_normal_nt0,head_binormal_nt0]; 
[xup1,xup2,xup3,xdown1,xdown2,xdown3] = Get_LF_Aligned_Head(B_nt0,neck_nt0);  % xup1: M1*M2.


    
    
%% Plot: 3D flow.

    fig = 8;
    figure(fig)
    clf; set(gcf,'color','w');
    hold on; axis equal; 
    view([125 20]); 
    %box off; axis off;
    xlabel('x'); ylabel('y'); zlabel('z')
    xlim([-1.53 xl_field(2)]); ylim([yl_field(1) 1.03]); zlim([-0.03 0.5]); 
    
    
    % Plot sperm.
    plot3(Xtail,Ytail,Ztail,'color',[36 100 171]./255,'linewidth',2);    %'color',[0.8500    0.3250    0.0980]
    CO_up(:,:,1) = 0.6*ones(size(xup1)); CO_up(:,:,2) = 0*ones(size(xup1)); CO_up(:,:,3) = 0*ones(size(xup1));
    CO_down(:,:,1) = 0*ones(size(xdown1)); CO_down(:,:,2) = 0.6*ones(size(xdown1)); CO_down(:,:,3) = 0*ones(size(xdown1));
    s1= surf(xup1,xup2,xup3,CO_up);
    s2= surf(xdown1,xdown2,xdown3,CO_down);
    s1.EdgeColor='none' ;
    s2.EdgeColor='none' ;
    
    % Plot 3D flow streamlines.
    [sx_3D,sy_3D,sz_3D] = meshgrid(xg(3:18),yg(6:18),zg(13)); 
    Generate_FlowStreamlines(Xg,Yg,Zg,Ug0,Vg0,Wg0,sx_3D,sy_3D,sz_3D, vel_dim);
        
    
    % Plot 2D slice planes.
    slice_height = 0:0.15:0.3;
    for i_slice = 1:length(slice_height)
        XL=[xl_field(1) xl_field(2)];
        YL=[yl_field(1) yl_field(2)];
        ZL=[slice_height(i_slice) zl_field(2)]; 
        [xgrid,ygrid,zgrid] = meshgrid(linspace(XL(1),XL(2),10),...
                    linspace(YL(1),YL(2),10),...
                    linspace(ZL(1),ZL(2),10)); 
        if i_slice == 1      
             CO(:,:,1) = ones(10)*.2;  
             CO(:,:,2) = ones(10)*.2; 
             CO(:,:,3) = ones(10)*0.2;  
        else
            CO(:,:,1) = ones(10)*.0; % red 
            CO(:,:,2) = ones(10)*.4470; % green 
            CO(:,:,3) = ones(10)*0.7410; % blue 
        end
        surfZ = surf(xgrid(:,:,1),ygrid(:,:,1),zgrid(:,:,1), CO,'FaceAlpha',0.1);  
        surfZ.EdgeColor='none';
    end
    shading interp
       
    camlight right
    lighting gouraud
    colorbar; caxis([0 20]);
  
    

   
%% Plot: 2D streamslices.


 for i_slice=1:length(slice_height)
    figure(fig+i_slice)
    clf;  set(gcf,'color','w');
    hold on; axis equal; view(2); 
    %box off; axis off;    
    xlabel('x'); ylabel('y'); zlabel('z')
    xlim([xl_field(1) xl_field(2)]); ylim([yl_field(1) yl_field(2)]); zlim([0 zl_field(2)]); 
    %xlim([-1.55 xl_field(2)]); ylim([yl_field(1) 1.05]); zlim([-0.03 0.5]); 
    
    % Plot sperm.
    plot3(Xtail,Ytail,Ztail,'color',[36 100 171]./255,'linewidth',2);    %'color',[0.8500    0.3250    0.0980]
    CO_up(:,:,1) = 0.6*ones(size(xup1)); CO_up(:,:,2) = 0*ones(size(xup1)); CO_up(:,:,3) = 0*ones(size(xup1));
    CO_down(:,:,1) = 0*ones(size(xdown1)); CO_down(:,:,2) = 0.6*ones(size(xdown1)); CO_down(:,:,3) = 0*ones(size(xdown1));
    s1= surf(xup1,xup2,xup3,CO_up);
    s2= surf(xdown1,xdown2,xdown3,CO_down);
    s1.EdgeColor='none' ;
    s2.EdgeColor='none' ;
    
   
    % Plot 2D streamslices.
    ss=streamslice(Xg,Yg,Zg,Ug0,Vg0,Wg0,[],[],[slice_height(i_slice)],0.3);
    set(ss,'color','white','linewidth',10);
    slice(Xg,Yg,Zg,vel_dim,[],[],[slice_height(i_slice)]);
    shading interp; 
    colormap(gca,'parula')
    colorbar; caxis([0 20]);
 end
  
end



%% For Fig 3.

if 1
    
fig = 1;
figure(fig)
    clf; set(gcf,'color','w');
    hold on; axis equal; 
    view([125 20]); 
    box off; axis off;
    xlim([-1.5 1]); ylim([-1 1]); zlim([-1 1]); 
    xlabel('x'); ylabel('y'); zlabel('z')
    
     % Plot Cartesian axes.
    O = [-0.5 -0.5 -0.5];
    Cx=[1; 0; 0]; 
    Cy=[0; 1; 0]; 
    Cz=[0; 0; 1]; 
    C=[Cx Cy Cz];
    Plot_Cartesian_Axes_1(fig,O,C)  
    
camlight right
lighting gouraud

end