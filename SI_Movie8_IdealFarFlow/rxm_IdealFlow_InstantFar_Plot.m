function rxm_IdealFlow_InstantFar_Plot

% This function plots the far-field average flow around the idealized sperm, relative to the comoving frame of reference.
% Figures:
% 1. 3D average flow. 
% 2. flagellar shape.
% 3. Colorbar.

clc


%% Fig. 1: 3D flow.

if 1
    %=====================================================================================================
    % Get flow field.
    load FarField_CF_alpha0.6.mat;
    %np_beat = 20;
    nt0 = 8; % selected instant moments: nt0 = [4 8 12 16 20].
       
    xg=linspace(xl_field(1),xl_field(2),Nx_field);
    yg=linspace(yl_field(1),yl_field(2),Ny_field); 
    zg=linspace(zl_field(1),zl_field(2),Nz_field); 
    [Xg,Yg,Zg]=meshgrid(xg,yg,zg);   
    
    Ug_nt0 = reshape(Ug(:,nt0),size(Xg,1),size(Xg,2),size(Xg,3));
    Vg_nt0 = reshape(Vg(:,nt0),size(Xg,1),size(Xg,2),size(Xg,3));
    Wg_nt0 = reshape(Wg(:,nt0),size(Xg,1),size(Xg,2),size(Xg,3));
    vel_raw=(Ug_nt0.*Ug_nt0+Vg_nt0.*Vg_nt0+Wg_nt0.*Wg_nt0).^0.5;
    vel_nd=vel_raw.*2*pi; % Transform the velocity units into 'flagellum length/ beat cycle'.

    
    %=====================================================================================================
    % Get sperm shape.
    load XNodes_Lab_alpha0.6.mat;
    %np_beat = 20;
    nt_XNodes=find(tRange==tRange_flow(end));
    nt0_XNodes=find(tRange==tRange_flow(nt0));
    Q=Nhh+Ns;   
    Xtail_raw = Xs(Nhh+1:Q,1:nt_XNodes);
    Ytail_raw = Xs(Q+Nhh+1:2*Q,1:nt_XNodes); 
    Ztail_raw = Xs(2*Q+Nhh+1:3*Q,1:nt_XNodes); 
     % Get the mean flagellar shape direction, based on which the instant sperm shape relative to the comoving frame can be obtained later by rotation.
    Xtail_0=Xtail_raw-repmat(Xtail_raw(1,:),Ns,1);
    Ytail_0=Ytail_raw-repmat(Ytail_raw(1,:),Ns,1);
    Ztail_0=Ztail_raw-repmat(Ztail_raw(1,:),Ns,1);   
    Xtail_mean= mean(Xtail_0,2); %Ns*1  %The mean flagellar shape relative to the lab frame.
    Ytail_mean= mean(Ytail_0,2);
    Ztail_mean= mean(Ztail_0,2);
    dir_new = [-1, 0, 0];
    dir_raw=[Xtail_mean(1)-Xtail_mean(end),Ytail_mean(1)-Ytail_mean(end),Ztail_mean(1)-Ztail_mean(end)];
    
    % Get the instant sperm shape.
    % Aligned tail.
    Xtail_nt0 = Xtail_0(:,nt0_XNodes); 
    Ytail_nt0 = Ytail_0(:,nt0_XNodes); 
    Ztail_nt0 = Ztail_0(:,nt0_XNodes); 
    [Xtail_align,Ytail_align,Ztail_align] = rxm_dir_align(dir_raw,dir_new,Xtail_nt0,Ytail_nt0,Ztail_nt0); %Xtail_align: Ns*1.
    % Aligned head.       
    head_tangent_nt0=head_tangent(nt0_XNodes,:)';
    head_normal_nt0=head_normal(nt0_XNodes,:)';
    head_binormal_nt0=cross(head_tangent_nt0,head_normal_nt0);
    B_nt0=[head_tangent_nt0,head_normal_nt0,head_binormal_nt0];
    [Bx_nt0_align,By_nt0_align,Bz_nt0_align] = rxm_dir_align(dir_raw,dir_new,B_nt0(1,:)',B_nt0(2,:)',B_nt0(3,:)');
    B_nt0_aligned = [Bx_nt0_align'; By_nt0_align'; Bz_nt0_align']; %3*3.  
    [xup1,xup2,xup3,xdown1,xdown2,xdown3,Ind_xup_neck] = Get_CF_Aligned_Head(B_nt0_aligned);  % xup1: M1*M2.   
    Neck =  [xup1(Ind_xup_neck) xup2(Ind_xup_neck) xup3(Ind_xup_neck)]; 
    Xtail = Xtail_align + Neck(1);  
    Ytail = Ytail_align + Neck(2);
    Ztail = Ztail_align + Neck(3);
    
     
    %================================================================================================================================================
    % Plot.
    fig = 1;
    
    figure(fig)
    clf; set(gcf,'color','w');
    hold on; axis equal; %axis off; 
    xlabel('x'); ylabel('y'); zlabel('z');
    %//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    XL=[xl_field(1)-2.5 xl_field(2)+2.5]; YL=[yl_field(1)-1.6 yl_field(2)]; ZL=[zl_field(1)-2.1 zl_field(2)+2.1];  % for cross-sectional flows.
    %XL=[0.5 1.2]; YL=[-0.5 0.2]; ZL=[-0.2 0.5];   % for XYZ axes. 
    xlim([XL(1) XL(2)]); ylim([YL(1) YL(2)]); zlim([ZL(1) ZL(2)]); 
    view([-35 20])   
    
if 1   
    % Plot sperm.
    plot3(Xtail,Ytail,Ztail,'color',[36 100 171]./255,'linewidth',2);    %'color',[0.8500    0.3250    0.0980]
    CO_up(:,:,1) = 0.6*ones(size(xup1)); CO_up(:,:,2) = 0*ones(size(xup1)); CO_up(:,:,3) = 0*ones(size(xup1));
    CO_down(:,:,1) = 0*ones(size(xdown1)); CO_down(:,:,2) = 0.6*ones(size(xdown1)); CO_down(:,:,3) = 0*ones(size(xdown1));
    s1= surf(xup1,xup2,xup3,CO_up);
    s2= surf(xdown1,xdown2,xdown3,CO_down);
    s1.EdgeColor='none' ;
    s2.EdgeColor='none' ;
    
    
    if 0 % Plot 2D slice planes.
    % X planes.
    XL1_mat = [-3 0 4];
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
    YL1_mat = [-5 0];
    for i_y = 1:length(YL1_mat)
        XL_temp=[XL(1) XL(2)];
        YL_temp=[YL1_mat(i_y) YL(2)];
        ZL_temp=[ZL(1) ZL(2)]; 
        [xgrid,ygrid,zgrid] = meshgrid(linspace(XL_temp(1),XL_temp(2),10),...
                    linspace(YL_temp(1),YL_temp(2),10),...
                    linspace(ZL_temp(1),ZL_temp(2),10));     
        CO(:,:,1) = ones(10)*.0; % red 
        CO(:,:,2) = ones(10)*.4470; % green 
        CO(:,:,3) = ones(10)*0.7410; % blue 
        surfY = surf(reshape(xgrid(1,:,:),[10 10]),reshape(ygrid(1,:,:),[10 10]),reshape(zgrid(1,:,:),[10,10]), CO,'FaceAlpha',0.1);  
        surfY.EdgeColor='none';
    end
    % Z plane.   
    ZL1_mat = [0 5];
    for i_z = 1:length(ZL1_mat)    
        XL_temp=[XL(1) XL(2)];
        YL_temp=[YL(1) YL(2)];
        ZL_temp=[ZL1_mat(i_z) ZL(2)]; 
        [xgrid,ygrid,zgrid] = meshgrid(linspace(XL_temp(1),XL_temp(2),10),...
                    linspace(YL_temp(1),YL_temp(2),10),...
                    linspace(ZL_temp(1),ZL_temp(2),10));     
        CO(:,:,1) = ones(10)*.0; % red 
        CO(:,:,2) = ones(10)*.4470; % green 
        CO(:,:,3) = ones(10)*0.7410; % blue 
        surfZ = surf(reshape(xgrid(:,:,1),[10 10]),reshape(ygrid(:,:,1),[10 10]),reshape(zgrid(:,:,1),[10,10]), CO,'FaceAlpha',0.1);  
        surfZ.EdgeColor='none';
    end
    end
    
  if 1
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////
    % Plot (cross-sectional) flow streamlines.
    % cross-sectional flow at x=0.
    [sx_3D,sy_3D,sz_3D] = meshgrid(xg(23),yg(1:end),zg(1:1:end));
    Generate_FlowStreamlines(Xg,Yg,Zg,Ug_nt0,Vg_nt0,Wg_nt0,sx_3D,sy_3D,sz_3D, vel_nd);
    % cross-sectional flow at y=0.
    [sx_3D,sy_3D,sz_3D] = meshgrid(xg(1:1:end),yg(21),zg(1:1:end));
    Generate_FlowStreamlines(Xg,Yg,Zg,Ug_nt0,Vg_nt0,Wg_nt0,sx_3D,sy_3D,sz_3D, vel_nd);
    % cross-sectional flow at z=0.
    [sx_3D,sy_3D,sz_3D] = meshgrid(xg(1:1:end),yg(1:1:end),zg(21));
    Generate_FlowStreamlines(Xg,Yg,Zg,Ug_nt0,Vg_nt0,Wg_nt0,sx_3D,sy_3D,sz_3D, vel_nd);
  end   
    shading interp
    colorbar; 
    caxis([0 0.0001]); 
        
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////
    % Arrows indicating the flow directions.
    if 1   
        % Note that the positions of the arrows need slight adjustments for
        % different moments. This is because, for one thing, some positions
        % would cause the arrows to be hidden by the dense streamlines,
        % while for the other, some clearly-seen positions may deviate too
        % much from the middle of the cross-sectional plane, leading to the
        % wrong pusher/ puller features.
        if nt0==0 %原版位置。
            [sx_cone1,sy_cone1,sz_cone1] = meshgrid(xg([23]),yg([3 19 40]),zg([4 18 40]));  % arrows located at the top.                
            [sx_cone2,sy_cone2,sz_cone2] = meshgrid(xg([23]),yg([3 19 40]),zg([2 18 41]));  % arrows located at the bottom.        
            [sx_cone3,sy_cone3,sz_cone3] = meshgrid(xg([10 23 40]),yg([3 26 39]),zg([22]));  % arrows located at the left.        
            [sx_cone4,sy_cone4,sz_cone4] = meshgrid(xg([10 23 40]),yg([3 26 39]),zg([22]));  % arrows located at the right.       
            [sx_cone5,sy_cone5,sz_cone5] = meshgrid(xg([2 23 36]),yg([19]),zg([8 22 40]));  % arrows located in the front.
            [sx_cone6,sy_cone6,sz_cone6] = meshgrid(xg([6 23 43]),yg([19]),zg([8 23 40]));  % arrows located at the back.
        else%if nt0==4  %标准中心位置。
            [sx_cone1,sy_cone1,sz_cone1] = meshgrid(xg([23]),yg([3 21 40]),zg([5 18 37]));  
            [sx_cone2,sy_cone2,sz_cone2] = meshgrid(xg([23]),yg([3 21 40]),zg([5 18 37]));  
            [sx_cone3,sy_cone3,sz_cone3] = meshgrid(xg([10 23 40]),yg([5 18 37]),zg([21]));  
            [sx_cone4,sy_cone4,sz_cone4] = meshgrid(xg([10 23 40]),yg([5 18 37]),zg([21]));  
            [sx_cone5,sy_cone5,sz_cone5] = meshgrid(xg([5 20 39]),yg([21]),zg([8 21 40]));  
            [sx_cone6,sy_cone6,sz_cone6] = meshgrid(xg([5 20 39]),yg([21]),zg([8 21 40]));          
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
        
        Ug_cone = Ug_nt0./vel_raw;  Vg_cone = Vg_nt0./vel_raw;  Wg_cone = Wg_nt0./vel_raw;
        scal = 8; %scal = 3; 
        Plot_StreamCone(fig,Xg,Yg,Zg,Ug_cone,Vg_cone,Wg_cone,sx_cone1,sy_cone1,sz_cone1,scal)
        Plot_StreamCone(fig,Xg,Yg,Zg,Ug_cone,Vg_cone,Wg_cone,sx_cone2,sy_cone2,sz_cone2,scal)
        Plot_StreamCone(fig,Xg,Yg,Zg,Ug_cone,Vg_cone,Wg_cone,sx_cone3,sy_cone3,sz_cone3,scal)
        Plot_StreamCone(fig,Xg,Yg,Zg,Ug_cone,Vg_cone,Wg_cone,sx_cone4,sy_cone4,sz_cone4,scal)
        Plot_StreamCone(fig,Xg,Yg,Zg,Ug_cone,Vg_cone,Wg_cone,sx_cone5,sy_cone5,sz_cone5,scal)
        Plot_StreamCone(fig,Xg,Yg,Zg,Ug_cone,Vg_cone,Wg_cone,sx_cone6,sy_cone6,sz_cone6,scal)        
   end
    
    camlight right
    lighting gouraud
    return

end


% Plot Cartesian axes.
    % Plot Cartesian axes.
    O = [1 0 0]; % for PCA modes.
    Cx=[-1; 0; 0]; 
    Cy=[0; -1; 0]; 
    Cz=[0; 0; 1]; 
    C=[Cx Cy Cz];
    Plot_Cartesian_Axes(fig,O,C)  
    camlight right
    lighting gouraud    
    return
end




%% Plot the flagellar sequence in its waveform.

if 0
    
%np_beat = 20;
nt0 = 20; % selected instant moments: nt0 = [4 8 12 16 20].
       
head_tangent = [-1; 0; 0];
head_normal = [0; -1; 0];
head_binormal = cross(head_tangent,head_normal);
B = [head_tangent,head_normal,head_binormal];
[xup1,xup2,xup3,xdown1,xdown2,xdown3,~] = Get_CF_Aligned_Head(B);  % xup1: M1*M2.

load IdealWF_alpha0.6.mat;
x = x + max(xup1(:)); %size(x)=[ns nt]

 
figure(3)
clf;hold on;set(gcf,'color','w')
view(-65,20); 
axis equal;
if 0
XL = [-0.1,1.2];YL = [-0.3,0.3];ZL = [-0.2,0.2]; 
xlim([XL(1),XL(2)]);ylim([YL(1),YL(2)]);zlim([ZL(1),ZL(2)]);
xticks([0 1]);yticks([-0.2 0 0.2]);zticks([-0.2 0 0.2]);  %xyz-model
xticklabels({'0','1'});yticklabels({'0.2','0','-0.2'});zticklabels({'-0.2','0','0.2'});  %xyz-model
ax = gca;ax.FontSize = 30; ax.TickLabelInterpreter = 'latex'; %ax.FontSize = 15; 
xlh=xlabel('$\xi_3$ ','interpreter','latex','FontSize',40);  %,'FontSize',20);
ylh=ylabel('$\xi_1$ ','interpreter','latex','FontSize',40);
zlh=zlabel('$\xi_2$ ','interpreter','latex','FontSize',40);
set(xlh,'position',[0.3 YL(1)-0.05 ZL(1)]);
set(ylh,'position',[XL(1)-0.05 0.1 ZL(1)]);
set(zlh,'position',[XL(1)-0.15 YL(2) 0.05]);
end
axis off; box off;

% Plot head
source_ligth = [90 80];
k = [0.9 1 1 5]; 
s1= surfl(xup1,xup2,xup3,source_ligth,k );
s2= surfl(xdown1,xdown2,xdown3,source_ligth,k );
s1.EdgeColor='none' ;
s1.FaceColor=  [0.6 0 0]; 
s2.EdgeColor='none' ;
s2.FaceColor=  [0 0.6 0];
camlight right
lighting gouraud    
    
% Plot tail waveform
for i_nt = [4 8 12 16 20]
   p = plot3(x(:,i_nt),y(:,i_nt),z(:,i_nt),'color',[36 100 171]./255,'linewidth',8);
   % transparency
   p.Color(4) = 0.2;
end
plot3(x(:,nt0),y(:,nt0),z(:,nt0),'color',[36 100 171]./255,'linewidth',8)


return

end



    



%% Fig. 3: Colorbar.

figure(3)
clf; set(gcf,'color','w');
hold on; axis equal; axis off; box off;
       
colormap(parula); 
clim = 0.0001;
caxis([0,clim]);

c = colorbar('Ticks',[0,clim],'FontSize',40,'TickLabels',{'0','$10^{-4}$'},'fontname','Times');  %'FontSize',65/ 80
c.Label.String = 'Flow velocity';
c.TickLabelInterpreter='latex';
c.Position = [0.7 0.1 0.005 0.4];  % orientation: vertical
return

c = colorbar('southoutside','Ticks',[0,clim],'TickLabels',{'0','$5\times10^{-5}$'},'FontSize',80,'FontName','Times');  
c.Label.String = 'Flow velocity';
c.TickLabelInterpreter='latex'; 
c.Position = [0.1 0.3 0.44 0.015];  % orientation: horizontal


