function rxm_IdealFlow_CF_Close_csYZ

% This function generates a video containing the instant flow streamlines
% of the idealized model's near flow.

clc

if 1
%================================================================================
% Plot the colorbar.
figure(1) 
clf; 
hold on; set(gcf,'color','w'); axis equal; %axis off;
     
fs = 30; 
caxis_max = 0.1;
    c=colorbar; 
    c.Position=[0.7 0.2 0.007 0.5];    
    caxis([0 caxis_max]); c.Ticks=[0 caxis_max]; 
    text(0.95,0.2,-0.18,{'Flow velocity'},'FontName','times','FontSize',fs,'rotation',90);  %'FontWeight','bold'  <-- for close-field flow video
    c.TickLabels={'0','0.1'};    
    c.FontName = 'times'; c.TickLabelInterpreter='latex'; c.FontSize=fs; 
 
return
end 
 


FieldRange = 'c'; %close


%% Data pre-processing.

% Get comoving-frame flow field.
if FieldRange == 'n'
    load NearFlow_CF_alpha0.6_Tfine.mat; 
elseif FieldRange == 'c'
    load CloseFlow_CF_alpha0.6_Tfine.mat; 
end
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
        

% Determine the time period for the video.    
nt_bc = 100;
nt_rev = 600; % This is because the '_Tfine' data have 100 time points for each beat cycle, and the case 'alpha0.6' is of the comoving revolution cycle with 6 beat cycles.

fs = 32;   
fig = 1;
figure(fig) 
clf; 
hold on;  

i_nt = 1*nt_bc;



%% Plot: comoving-frame sperm body snapshots.

if 0
  for i_nt_temp = [nt_bc:nt_bc:nt_rev]
    % Get the instant sperm shape.
    nt0_X = find(tRange==tRange_flow(i_nt_temp));
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
    
   
    set(gcf,'color','w'); axis equal; box on; box off; axis off;
    view([-35 15]);  %view([-55 15]);
    XL = [-0.05 0.95];  YL = [-0.14 0.14];  ZL = [-0.14 0.14]; 
    xlim([XL(1) XL(2)]); ylim([YL(1) YL(2)]); zlim([ZL(1) ZL(2)]);  
    xlabel('x'); ylabel('y'); zlabel('z');
    %================================================================================
    % Plot sperm. 
    p_tail = plot3(Xtail_CF,Ytail_CF,Ztail_CF,'color',[36 100 171]./255,'linewidth',5);    %'color',[0.8500    0.3250    0.0980]
    CO_up(:,:,1) = 0.6*ones(size(xup1)); CO_up(:,:,2) = 0*ones(size(xup1)); CO_up(:,:,3) = 0*ones(size(xup1));
    CO_down(:,:,1) = 0*ones(size(xdown1)); CO_down(:,:,2) = 0.6*ones(size(xdown1)); CO_down(:,:,3) = 0*ones(size(xdown1));
    s1= surf(xup1,xup2,xup3,CO_up);
    s2= surf(xdown1,xdown2,xdown3,CO_down);
    s1.EdgeColor='none' ;
    s2.EdgeColor='none' ;
   if  i_nt_temp==i_nt
    s1.FaceAlpha = 1;  
    s2.FaceAlpha = 1;
    p_tail.Color(4) = 1;
   else
    s1.FaceAlpha = 0.1;  
    s2.FaceAlpha = 0.1; 
    p_tail.Color(4) = 0.2;
   end
  end
    %================================================================================
    % Plot cross-sectional planes: YZ-planes. 
        XL1_mat = 0.9; %[0 0.5 0.9];   %[0 0.52 0.93];  %
        for i_x = 1:length(XL1_mat)
            XL_temp=[XL1_mat(i_x) XL(2)];
            YL_temp=[YL(1) YL(2)];
            ZL_temp=[ZL(1) ZL(2)]; 
            [xgrid,ygrid,zgrid] = meshgrid(linspace(XL_temp(1),XL_temp(2),10),...
                    linspace(YL_temp(1),YL_temp(2),10),...
                    linspace(ZL_temp(1),ZL_temp(2),10));     
            if 1 % plane color: grey-blue
                CO(:,:,1) = ones(10)*.0; % red 
                CO(:,:,2) = ones(10)*.4470; % green 
                CO(:,:,3) = ones(10)*0.7410; % blue 
            else % plane color: brown-red
                CO(:,:,1) = ones(10)*.85; % red
                CO(:,:,2) = ones(10)*0.3250 ; % green
                CO(:,:,3) = ones(10)*0.0980; % blue
            end
            surfX = surf(reshape(xgrid(:,1,:),[10 10]),reshape(ygrid(:,1,:),[10 10]),reshape(zgrid(:,1,:),[10,10]), CO,'FaceAlpha',0.1);   %,'FaceAlpha',0.1
            surfX.EdgeColor='none';
        end
    camlight right
    lighting gouraud  
    
    
elseif 1
%% Plot: comoving-frame flow.

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
    
    
%Get instant flow.
Ug_nt0 = reshape(Ug(:,i_nt),size(Xg,1),size(Xg,2),size(Xg,3));
Vg_nt0 = reshape(Vg(:,i_nt),size(Xg,1),size(Xg,2),size(Xg,3));
Wg_nt0 = reshape(Wg(:,i_nt),size(Xg,1),size(Xg,2),size(Xg,3));
vel_raw=(Ug_nt0.*Ug_nt0+Vg_nt0.*Vg_nt0+Wg_nt0.*Wg_nt0).^0.5;
vel_nd=vel_raw.*2*pi; % Transform the velocity units into 'flagellum length/ beat cycle'.    
    
if FieldRange == 'n'
    caxis_max = 0.01; %<-- for near-field flow video
elseif FieldRange == 'c'
    caxis_max =  0.1; %<-- for close-field flow video
end

if 0
%================================================================================
%% Plot: Instantaneous flow. cs-YZ-plane streamlines -- near head (x=0).  end flagellum (x=0.9).  
%================================================================================
    set(gcf,'color','w'); axis equal; view([-90 0]);  
    xlim([xl_field(1) xl_field(2)]); ylim([yl_field(1) yl_field(2)]); zlim([zl_field(1) zl_field(2)]); 
    xticks([]); yticks([]); zticks([]);
    %ylh=ylabel('$y_c$','interpreter','latex','FontSize',1.5*fs);
    %zlh=zlabel('$z_c$','interpreter','latex','FontSize',1.5*fs,'Rotation',0);
    %set(ylh,'position',[0,0,-1]);  set(zlh,'position',[0,1.3,0])
%================================================================================
% Plot streamslice.
    ss=streamslice(Xg,Yg,Zg,Ug_nt0,Vg_nt0,Wg_nt0,[0.9],[],[],0.3);
    set(ss,'color','white','linewidth',12);   
    slice(Xg,Yg,Zg,vel_nd,[0.9],[],[]);
    shading interp; 
    colormap(gca,'parula'); caxis([0 caxis_max]);
    % Plot sperm.
    plot3(Xtail_CF,Ytail_CF,Ztail_CF,'color',[36 100 171]./255,'linewidth',5);    %'color',[0.8500    0.3250    0.0980]
    CO_up(:,:,1) = 0.6*ones(size(xup1)); CO_up(:,:,2) = 0*ones(size(xup1)); CO_up(:,:,3) = 0*ones(size(xup1));
    CO_down(:,:,1) = 0*ones(size(xdown1)); CO_down(:,:,2) = 0.6*ones(size(xdown1)); CO_down(:,:,3) = 0*ones(size(xdown1));
    s1= surf(xup1,xup2,xup3,CO_up);
    s2= surf(xdown1,xdown2,xdown3,CO_down);
    s1.EdgeColor='none' ;
    s2.EdgeColor='none' ;
    
elseif 0           
%================================================================================
%% Plot: Cumulative time-averaged flow. cs-YZ-plane streamlines -- near head (x=0).  end flagellum (x=0.9).  
%================================================================================

%Get time-averaged flow.
Ug_mean = mean(Ug(:,1:i_nt),2);  
Vg_mean = mean(Vg(:,1:i_nt),2);   
Wg_mean = mean(Wg(:,1:i_nt),2);   
Ug_mean = reshape(Ug_mean,size(Xg,1),size(Xg,2),size(Xg,3));
Vg_mean = reshape(Vg_mean,size(Xg,1),size(Xg,2),size(Xg,3));
Wg_mean = reshape(Wg_mean,size(Xg,1),size(Xg,2),size(Xg,3)); 
vel_mean_raw=(Ug_mean.*Ug_mean+Vg_mean.*Vg_mean+Wg_mean.*Wg_mean).^0.5;
vel_mean_nd=vel_mean_raw.*2*pi; % Transform the velocity units into 'flagellum length/ beat cycle'.    
    
%================================================================================
    set(gcf,'color','w'); axis equal; view([-90 0]);  
    xlim([xl_field(1) xl_field(2)]); ylim([yl_field(1) yl_field(2)]); zlim([zl_field(1) zl_field(2)]); 
    xticks([]); yticks([]); zticks([]);
    %ylh=ylabel('$y_c$','interpreter','latex','FontSize',1.5*fs);
    %zlh=zlabel('$z_c$','interpreter','latex','FontSize',1.5*fs,'Rotation',0);
    %set(ylh,'position',[0,0,-1]);  set(zlh,'position',[0,1.3,0])
%================================================================================
% Plot streamslice.
    ss=streamslice(Xg,Yg,Zg,Ug_mean,Vg_mean,Wg_mean,[0.9],[],[],0.3);
    set(ss,'color','white','linewidth',12);   
    slice(Xg,Yg,Zg,vel_mean_nd,[0.9],[],[]);
    shading interp; 
    colormap(gca,'parula'); caxis([0 caxis_max]);
    % Plot sperm.
    plot3(Xtail_CF,Ytail_CF,Ztail_CF,'color',[36 100 171]./255,'linewidth',5);    %'color',[0.8500    0.3250    0.0980]
    CO_up(:,:,1) = 0.6*ones(size(xup1)); CO_up(:,:,2) = 0*ones(size(xup1)); CO_up(:,:,3) = 0*ones(size(xup1));
    CO_down(:,:,1) = 0*ones(size(xdown1)); CO_down(:,:,2) = 0.6*ones(size(xdown1)); CO_down(:,:,3) = 0*ones(size(xdown1));
    s1= surf(xup1,xup2,xup3,CO_up);
    s2= surf(xdown1,xdown2,xdown3,CO_down);
    s1.EdgeColor='none' ;
    s2.EdgeColor='none' ;

end      
end


