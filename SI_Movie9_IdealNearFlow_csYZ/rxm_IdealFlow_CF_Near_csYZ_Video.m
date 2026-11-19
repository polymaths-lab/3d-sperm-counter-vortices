function rxm_IdealFlow_CF_Near_csYZ_Video

% This function generates a video containing the instant flow streamlines
% of the idealized model's near flow.

clc


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
        

%% Generate the video. 

% Determine the time period for the video.    
nt_bc = 100;
nt_rev = 600; % This is because the '_Tfine' data have 100 time points for each beat cycle, and the case 'alpha0.6' is of the comoving revolution cycle with 6 beat cycles.


fs = 32;   
fig = 2;
figure(fig) 
if FieldRange == 'n'
    v=VideoWriter('NearFlow_alpha0.6_csYZ_test1','MPEG-4');
elseif FieldRange == 'c'
    v=VideoWriter('CloseFlow_alpha0.6_csYZ_test1','MPEG-4');
end
v.FrameRate=20;  
open(v);     
for i_nt = 1:nt_rev+1
    i_nt
    
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
    
    
clf; 
hold on;  

%%
%================================================================================
% First row: comoving-frame sperm body motion + flagellar point tracks.
%================================================================================
    set(gca,'position',[0.35 0.65 0.6 0.35]);  %set(gca,'position',[0.25 0.5 0.6 0.5]);  %get(gca,'position')    
    set(gcf,'color','w'); axis equal; box on; box off; axis off;
    view([-30 15]);  %view([-45 20]);
    XL = [-0.05 0.95];  YL = [-0.14 0.14];  ZL = [-0.14 0.14]; 
    xlim([XL(1) XL(2)]); ylim([YL(1) YL(2)]); zlim([ZL(1) ZL(2)]);  
    xlabel('x'); ylabel('y'); zlabel('z');
    %================================================================================
    % Plot sperm. 
    plot3(Xtail_CF,Ytail_CF,Ztail_CF,'color',[36 100 171]./255,'linewidth',2);    %'color',[0.8500    0.3250    0.0980]
    CO_up(:,:,1) = 0.6*ones(size(xup1)); CO_up(:,:,2) = 0*ones(size(xup1)); CO_up(:,:,3) = 0*ones(size(xup1));
    CO_down(:,:,1) = 0*ones(size(xdown1)); CO_down(:,:,2) = 0.6*ones(size(xdown1)); CO_down(:,:,3) = 0*ones(size(xdown1));
    s1= surf(xup1,xup2,xup3,CO_up);
    s2= surf(xdown1,xdown2,xdown3,CO_down);
    s1.EdgeColor='none' ;
    s2.EdgeColor='none' ;
    %================================================================================
    % Plot flagellar point traces. 
    pp1 = [Xtail_CF(49); Ytail_CF(49); Ztail_CF(49)];  %[Xtail_CF(Ns/2); Ytail_CF(Ns/2); Ztail_CF(Ns/2)]; %point position
    pp2 = [Xtail_CF(97); Ytail_CF(97); Ztail_CF(97)];  %[Xtail_CF(Ns); Ytail_CF(Ns); Ztail_CF(Ns)];
    scatter3(pp1(1),pp1(2),pp1(3),150,'k','filled');
    scatter3(pp2(1),pp2(2),pp2(3),150,'k','filled');
    if i_nt==1
        pt1 = pp1;
        pt2 = pp2;
    else
        pt1 = [pt1 pp1]; % 3*i_nt
        pt2 = [pt2 pp2]; 
        if 0 % trace plot is color-encoded
            u1 = Function_derivative_of_time(pt1,tRange(1:i_nt)./(2*pi),'time'); % 3*i_nt
            u2 = Function_derivative_of_time(pt2,tRange(1:i_nt)./(2*pi),'time');        
            u1_mag = sqrt(sum(u1.^2, 1)); % 1*i_nt
            u2_mag = sqrt(sum(u2.^2, 1));
            %
            pt1_plot = pt1;
            pt2_plot = pt2;
            pt1_plot(2,end) = NaN;
            pt2_plot(2,end) = NaN;
            patch(pt1_plot(1,:),pt1_plot(2,:),pt1_plot(3,:),u1_mag,'EdgeColor','interp','linewidth',3) 
            patch(pt2_plot(1,:),pt2_plot(2,:),pt2_plot(3,:),u2_mag,'EdgeColor','interp','linewidth',3) 
            colormap(parula)
        else % trace plot no varying color.
            plot3(pt1(1,:),pt1(2,:),pt1(3,:),'k','linewidth',3);
            plot3(pt2(1,:),pt2(2,:),pt2(3,:),'k','linewidth',3);
        end
    end
    %================================================================================
    % Plot cross-sectional planes: YZ-planes. 
        XL1_mat = [0 0.5 0.9];   %[0 0.52 0.93];  %
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
    %Texts indicating the three cross sections.
    text(-0.04,0.12,0.13,{'I'},'FontName','times','interpreter','latex','FontSize',fs);  %yz-plane: near head (x=0)  
    text(0.44,0.12,0.13,{'II'},'FontName','times','interpreter','latex','FontSize',fs);  %yz-plane: mid flagellum (x=0.5)  
    text(0.82,0.12,0.11,{'III'},'FontName','times','interpreter','latex','FontSize',fs);  %yz-plane: distal flagellum (x=0.9)  
    %================================================================================     
    % Plot the Cartesian axes for the comoving frame.  
    axes('Position',[0.3 0.75  0.1 0.1]);     %axes('Position',[0.13 0.6  0.1 0.1]);     
    set(gcf,'color','w');axis equal; box on; box off; axis off;
    view([-30 15]);  %view([-45 20]);
    xlabel('x'); ylabel('y'); zlabel('z');
    O =[0.1 -0.3 -0.25];  Cx=[-1; 0; 0]; Cy=[0; -1; 0]; Cz=[0; 0; 1]; C=[Cx Cy Cz];
    Plot_Cartesian_Axes_V2(fig,O,C) 
    % Text: axes labels.
    text(-1.85,-0.3,-0.25,{'$x_c$'},'FontName','times','interpreter','latex','FontSize',fs,'FontWeight','bold');  
    text(0.1,-1.7,-0.25,{'$y_c$'},'FontName','times','interpreter','latex','FontSize',fs,'FontWeight','bold');   
    text(0,0,1.3,{'$z_c$'},'FontName','times','interpreter','latex','FontSize',fs,'FontWeight','bold');
    % Text: sub-title.
    %text(-1.5,0,3.6,{'Comoving frame beat'},'FontName','times','FontSize',fs,'FontWeight','bold');   
    camlight right
    lighting gouraud  
    %================================================================================       
    % Text: varying time moment. 
    axes('Position',[0.01  0.8  0.1 0.1]);     
    set(gcf,'color','w');axis equal; box on; box off; axis off;
    t_nd_flow_nt0 = tRange_flow(i_nt)/(2*pi);
    instant_nt0 = roundn(t_nd_flow_nt0,-3);
    text(0.9,0.3,1.3,{['Time: ',num2str(instant_nt0)]},'FontName','times','FontSize',fs+10,'FontWeight','bold');  % time is in units of 'beat cycles'.  

    
    
%% Second row: instantaneous flow.

%Get instant flow.
Ug_nt0 = reshape(Ug(:,i_nt),size(Xg,1),size(Xg,2),size(Xg,3));
Vg_nt0 = reshape(Vg(:,i_nt),size(Xg,1),size(Xg,2),size(Xg,3));
Wg_nt0 = reshape(Wg(:,i_nt),size(Xg,1),size(Xg,2),size(Xg,3));
vel_raw=(Ug_nt0.*Ug_nt0+Vg_nt0.*Vg_nt0+Wg_nt0.*Wg_nt0).^0.5;
vel_nd=vel_raw.*2*pi; % Transform the velocity units into 'flagellum length/ beat cycle'.    
    
ps_2d_x = 0.2; %0.01;  
ps_2d_y = 0.37; %0.1;  
ps_2d_w = 0.25;  %0.35; 
ps_2d_h = 0.25;  %0.35;  
ps_2d_sp = -0.02; % 0.06; 

if FieldRange == 'n'
    caxis_max = 0.01; %<-- for near-field flow video
elseif FieldRange == 'c'
    caxis_max =  0.1; %<-- for close-field flow video
end
%================================================================================
% cs-YZ-plane streamlines -- near head (x=0).  
%================================================================================
    axes('Position',[ps_2d_x ps_2d_y ps_2d_w ps_2d_h]);
    hold on; set(gcf,'color','w'); axis equal; view([-90 0]);  
    xlim([xl_field(1) xl_field(2)]); ylim([yl_field(1) yl_field(2)]); zlim([zl_field(1) zl_field(2)]); 
    xticks([]); yticks([]); zticks([]);
    %ylh=ylabel('$y_c$','interpreter','latex','FontSize',1.5*fs);
    %zlh=zlabel('$z_c$','interpreter','latex','FontSize',1.5*fs,'Rotation',0);
    %set(ylh,'position',[0,0,-1]);  set(zlh,'position',[0,1.3,0])
    if FieldRange == 'n'
        text(0,1.2,1,{'I'},'FontName','times','interpreter','latex','FontSize',fs);  
    elseif FieldRange == 'c'
        text(0,0.24,0.18,{'I'},'FontName','times','interpreter','latex','FontSize',fs);  
    end
%================================================================================
% Plot streamslice.
    ss=streamslice(Xg,Yg,Zg,Ug_nt0,Vg_nt0,Wg_nt0,[0],[],[],0.5);
    set(ss,'color','white','linewidth',3);   
    slice(Xg,Yg,Zg,vel_nd,[0],[],[]);
    shading interp; 
    colormap(gca,'parula'); caxis([0 caxis_max]);
    % Plot sperm.
    plot3(Xtail_CF,Ytail_CF,Ztail_CF,'color',[36 100 171]./255,'linewidth',2);    %'color',[0.8500    0.3250    0.0980]
    CO_up(:,:,1) = 0.6*ones(size(xup1)); CO_up(:,:,2) = 0*ones(size(xup1)); CO_up(:,:,3) = 0*ones(size(xup1));
    CO_down(:,:,1) = 0*ones(size(xdown1)); CO_down(:,:,2) = 0.6*ones(size(xdown1)); CO_down(:,:,3) = 0*ones(size(xdown1));
    s1= surf(xup1,xup2,xup3,CO_up);
    s2= surf(xdown1,xdown2,xdown3,CO_down);
    s1.EdgeColor='none' ;
    s2.EdgeColor='none' ;
    
    

%================================================================================
% cs-YZ-plane streamlines -- middle flagellum (x=0.5).  
%================================================================================
    axes('Position',[ps_2d_x+ps_2d_w+ps_2d_sp ps_2d_y ps_2d_w ps_2d_h]);
    hold on; set(gcf,'color','w'); axis equal; view([-90 0]);  
    xlim([xl_field(1) xl_field(2)]); ylim([yl_field(1) yl_field(2)]); zlim([zl_field(1) zl_field(2)]); 
    xticks([]); yticks([]); zticks([]);
    %ylh=ylabel('$y_c$','interpreter','latex','FontSize',1.5*fs);
    %zlh=zlabel('$z_c$','interpreter','latex','FontSize',1.5*fs,'Rotation',0);
    %set(ylh,'position',[0,0,-1]);  set(zlh,'position',[0,1.3,0])
    if FieldRange == 'n'
        text(0,1.2,1,{'II'},'FontName','times','interpreter','latex','FontSize',fs);  
    elseif FieldRange == 'c'
        text(0,0.26,0.18,{'II'},'FontName','times','interpreter','latex','FontSize',fs);  
    end
%================================================================================
% Plot streamslice.
    ss=streamslice(Xg,Yg,Zg,Ug_nt0,Vg_nt0,Wg_nt0,[0.5],[],[],0.5);
    set(ss,'color','white','linewidth',3);   
    slice(Xg,Yg,Zg,vel_nd,[0.5],[],[]);
    shading interp; 
    colormap(gca,'parula'); caxis([0 caxis_max]);
    % Plot sperm.
    plot3(Xtail_CF,Ytail_CF,Ztail_CF,'color',[36 100 171]./255,'linewidth',2);    %'color',[0.8500    0.3250    0.0980]
    CO_up(:,:,1) = 0.6*ones(size(xup1)); CO_up(:,:,2) = 0*ones(size(xup1)); CO_up(:,:,3) = 0*ones(size(xup1));
    CO_down(:,:,1) = 0*ones(size(xdown1)); CO_down(:,:,2) = 0.6*ones(size(xdown1)); CO_down(:,:,3) = 0*ones(size(xdown1));
    s1= surf(xup1,xup2,xup3,CO_up);
    s2= surf(xdown1,xdown2,xdown3,CO_down);
    s1.EdgeColor='none' ;
    s2.EdgeColor='none' ;

    
%================================================================================
% cs-YZ-plane streamlines -- end flagellum (x=0.9).  
%================================================================================
    axes('Position',[ps_2d_x+2*ps_2d_w+2*ps_2d_sp ps_2d_y ps_2d_w ps_2d_h]);
    hold on; set(gcf,'color','w'); axis equal; view([-90 0]);  
    xlim([xl_field(1) xl_field(2)]); ylim([yl_field(1) yl_field(2)]); zlim([zl_field(1) zl_field(2)]); 
    xticks([]); yticks([]); zticks([]);
    %ylh=ylabel('$y_c$','interpreter','latex','FontSize',1.5*fs);
    %zlh=zlabel('$z_c$','interpreter','latex','FontSize',1.5*fs,'Rotation',0);
    %set(ylh,'position',[0,0,-1]);  set(zlh,'position',[0,1.3,0])
    if FieldRange == 'n'
        text(0,1.4,1,{'III'},'FontName','times','interpreter','latex','FontSize',fs);  
    elseif FieldRange == 'c'
        text(0,0.28,0.18,{'III'},'FontName','times','interpreter','latex','FontSize',fs);  
    end
%================================================================================
% Plot streamslice.
    ss=streamslice(Xg,Yg,Zg,Ug_nt0,Vg_nt0,Wg_nt0,[0.9],[],[],0.5);
    set(ss,'color','white','linewidth',3);   
    slice(Xg,Yg,Zg,vel_nd,[0.9],[],[]);
    shading interp; 
    colormap(gca,'parula'); caxis([0 caxis_max]);
    % Plot sperm.
    plot3(Xtail_CF,Ytail_CF,Ztail_CF,'color',[36 100 171]./255,'linewidth',2);    %'color',[0.8500    0.3250    0.0980]
    CO_up(:,:,1) = 0.6*ones(size(xup1)); CO_up(:,:,2) = 0*ones(size(xup1)); CO_up(:,:,3) = 0*ones(size(xup1));
    CO_down(:,:,1) = 0*ones(size(xdown1)); CO_down(:,:,2) = 0.6*ones(size(xdown1)); CO_down(:,:,3) = 0*ones(size(xdown1));
    s1= surf(xup1,xup2,xup3,CO_up);
    s2= surf(xdown1,xdown2,xdown3,CO_down);
    s1.EdgeColor='none' ;
    s2.EdgeColor='none' ;
%================================================================================
% Plot the colorbar.
    c=colorbar; 
    c.Position=[ps_2d_x+3*ps_2d_w+3*ps_2d_sp ps_2d_y 0.007 ps_2d_h];    
    caxis([0 caxis_max]); c.Ticks=[0 caxis_max]; 
    if FieldRange == 'n'
        text(1,-2,-0.9,{'Flow velocity'},'FontName','times','FontSize',fs,'rotation',90); %'FontWeight','bold' <-- for near-field flow video
        c.TickLabels={'0','0.01'};  
    elseif FieldRange == 'c'
        text(1,-0.4,-0.18,{'Flow velocity'},'FontName','times','FontSize',fs,'rotation',90);  %'FontWeight','bold'  <-- for close-field flow video
        c.TickLabels={'0','0.1'};    
    end
    c.FontName = 'times'; c.TickLabelInterpreter='latex'; c.FontSize=fs; 
  
 
%================================================================================       
% Text: instantaneous flow. 
    axes('Position',[0.0  0.5  0.1 0.1]);     
    set(gcf,'color','w');axis equal; box on; box off; axis off;
    text(0,0,0,{['Instantaneous flow']},'FontName','times','FontSize',fs,'FontWeight','bold');   
    
    
            
%% Third row: Cumulative time-averaged flow.

%Get time-averaged flow.
Ug_mean = mean(Ug(:,1:i_nt),2);  
Vg_mean = mean(Vg(:,1:i_nt),2);   
Wg_mean = mean(Wg(:,1:i_nt),2);   
Ug_mean = reshape(Ug_mean,size(Xg,1),size(Xg,2),size(Xg,3));
Vg_mean = reshape(Vg_mean,size(Xg,1),size(Xg,2),size(Xg,3));
Wg_mean = reshape(Wg_mean,size(Xg,1),size(Xg,2),size(Xg,3)); 
vel_mean_raw=(Ug_mean.*Ug_mean+Vg_mean.*Vg_mean+Wg_mean.*Wg_mean).^0.5;
vel_mean_nd=vel_mean_raw.*2*pi; % Transform the velocity units into 'flagellum length/ beat cycle'.    
    
ps_2d_x = 0.2;  
ps_2d_y = 0.05; 
ps_2d_w = 0.25;   
ps_2d_h = 0.25;    
ps_2d_sp = -0.02; 

%================================================================================
% cs-YZ-plane streamlines -- near head (x=0).  
%================================================================================
    axes('Position',[ps_2d_x ps_2d_y ps_2d_w ps_2d_h]);
    hold on; set(gcf,'color','w'); axis equal; view([-90 0]);  
    xlim([xl_field(1) xl_field(2)]); ylim([yl_field(1) yl_field(2)]); zlim([zl_field(1) zl_field(2)]); 
    xticks([]); yticks([]); zticks([]);
    %ylh=ylabel('$y_c$','interpreter','latex','FontSize',1.5*fs);
    %zlh=zlabel('$z_c$','interpreter','latex','FontSize',1.5*fs,'Rotation',0);
    %set(ylh,'position',[0,0,-1]);  set(zlh,'position',[0,1.3,0])
    if FieldRange == 'n'
        text(0,1.2,1,{'I'},'FontName','times','interpreter','latex','FontSize',fs);  
    elseif FieldRange == 'c'
        text(0,0.24,0.18,{'I'},'FontName','times','interpreter','latex','FontSize',fs);  
    end
%================================================================================
% Plot streamslice.
    ss=streamslice(Xg,Yg,Zg,Ug_mean,Vg_mean,Wg_mean,[0],[],[],0.5);
    set(ss,'color','white','linewidth',3);   
    slice(Xg,Yg,Zg,vel_mean_nd,[0],[],[]);
    shading interp; 
    colormap(gca,'parula'); caxis([0 caxis_max]);
    % Plot sperm.
    plot3(Xtail_CF,Ytail_CF,Ztail_CF,'color',[36 100 171]./255,'linewidth',2);    %'color',[0.8500    0.3250    0.0980]
    CO_up(:,:,1) = 0.6*ones(size(xup1)); CO_up(:,:,2) = 0*ones(size(xup1)); CO_up(:,:,3) = 0*ones(size(xup1));
    CO_down(:,:,1) = 0*ones(size(xdown1)); CO_down(:,:,2) = 0.6*ones(size(xdown1)); CO_down(:,:,3) = 0*ones(size(xdown1));
    s1= surf(xup1,xup2,xup3,CO_up);
    s2= surf(xdown1,xdown2,xdown3,CO_down);
    s1.EdgeColor='none' ;
    s2.EdgeColor='none' ;
    
    
%================================================================================
% cs-YZ-plane streamlines -- middle flagellum (x=0.5).  
%================================================================================
    axes('Position',[ps_2d_x+ps_2d_w+ps_2d_sp ps_2d_y ps_2d_w ps_2d_h]);
    hold on; set(gcf,'color','w'); axis equal; view([-90 0]);  
    xlim([xl_field(1) xl_field(2)]); ylim([yl_field(1) yl_field(2)]); zlim([zl_field(1) zl_field(2)]); 
    xticks([]); yticks([]); zticks([]);
    %ylh=ylabel('$y_c$','interpreter','latex','FontSize',1.5*fs);
    %zlh=zlabel('$z_c$','interpreter','latex','FontSize',1.5*fs,'Rotation',0);
    %set(ylh,'position',[0,0,-1]);  set(zlh,'position',[0,1.3,0])
    if FieldRange == 'n'
        text(0,1.3,1,{'II'},'FontName','times','interpreter','latex','FontSize',fs);  
    elseif FieldRange == 'c'
        text(0,0.26,0.18,{'II'},'FontName','times','interpreter','latex','FontSize',fs);  
    end
%================================================================================
% Plot streamslice.
    ss=streamslice(Xg,Yg,Zg,Ug_mean,Vg_mean,Wg_mean,[0.5],[],[],0.5);
    set(ss,'color','white','linewidth',3);   
    slice(Xg,Yg,Zg,vel_mean_nd,[0.5],[],[]);
    shading interp; 
    colormap(gca,'parula'); caxis([0 caxis_max]);
    % Plot sperm.
    plot3(Xtail_CF,Ytail_CF,Ztail_CF,'color',[36 100 171]./255,'linewidth',2);    %'color',[0.8500    0.3250    0.0980]
    CO_up(:,:,1) = 0.6*ones(size(xup1)); CO_up(:,:,2) = 0*ones(size(xup1)); CO_up(:,:,3) = 0*ones(size(xup1));
    CO_down(:,:,1) = 0*ones(size(xdown1)); CO_down(:,:,2) = 0.6*ones(size(xdown1)); CO_down(:,:,3) = 0*ones(size(xdown1));
    s1= surf(xup1,xup2,xup3,CO_up);
    s2= surf(xdown1,xdown2,xdown3,CO_down);
    s1.EdgeColor='none' ;
    s2.EdgeColor='none' ;
 
    
%================================================================================
% cs-YZ-plane streamlines -- end flagellum (x=0.9).  
%================================================================================
    axes('Position',[ps_2d_x+2*ps_2d_w+2*ps_2d_sp ps_2d_y ps_2d_w ps_2d_h]);
    hold on; set(gcf,'color','w'); axis equal; view([-90 0]);  
    xlim([xl_field(1) xl_field(2)]); ylim([yl_field(1) yl_field(2)]); zlim([zl_field(1) zl_field(2)]); 
    xticks([]); yticks([]); zticks([]);
    %ylh=ylabel('$y_c$','interpreter','latex','FontSize',1.5*fs);
    %zlh=zlabel('$z_c$','interpreter','latex','FontSize',1.5*fs,'Rotation',0);
    %set(ylh,'position',[0,0,-1]);  set(zlh,'position',[0,1.3,0])
    if FieldRange == 'n'
        text(0,1.4,1,{'III'},'FontName','times','interpreter','latex','FontSize',fs);  
    elseif FieldRange == 'c'
        text(0,0.28,0.18,{'III'},'FontName','times','interpreter','latex','FontSize',fs);  
    end
%================================================================================
% Plot streamslice.
    ss=streamslice(Xg,Yg,Zg,Ug_mean,Vg_mean,Wg_mean,[0.9],[],[],0.5);
    set(ss,'color','white','linewidth',3);   
    slice(Xg,Yg,Zg,vel_mean_nd,[0.9],[],[]);
    shading interp; 
    colormap(gca,'parula'); caxis([0 caxis_max]);
    % Plot sperm.
    plot3(Xtail_CF,Ytail_CF,Ztail_CF,'color',[36 100 171]./255,'linewidth',2);    %'color',[0.8500    0.3250    0.0980]
    CO_up(:,:,1) = 0.6*ones(size(xup1)); CO_up(:,:,2) = 0*ones(size(xup1)); CO_up(:,:,3) = 0*ones(size(xup1));
    CO_down(:,:,1) = 0*ones(size(xdown1)); CO_down(:,:,2) = 0.6*ones(size(xdown1)); CO_down(:,:,3) = 0*ones(size(xdown1));
    s1= surf(xup1,xup2,xup3,CO_up);
    s2= surf(xdown1,xdown2,xdown3,CO_down);
    s1.EdgeColor='none' ;
    s2.EdgeColor='none' ;
%================================================================================
% Plot the colorbar.
    c=colorbar; 
    c.Position=[ps_2d_x+3*ps_2d_w+3*ps_2d_sp ps_2d_y 0.007 ps_2d_h];    
    caxis([0 caxis_max]);    c.Ticks=[0 caxis_max]; 
    if FieldRange == 'n'
        text(1,-2,-0.9,{'Flow velocity'},'FontName','times','FontSize',fs,'rotation',90); %'FontWeight','bold' <-- for near-field flow video
        c.TickLabels={'0','0.01'};  
    elseif FieldRange == 'c'
        text(1,-0.4,-0.18,{'Flow velocity'},'FontName','times','FontSize',fs,'rotation',90);  %'FontWeight','bold'  <-- for close-field flow video
        c.TickLabels={'0','0.1'};    
    end
    c.FontName = 'times'; c.TickLabelInterpreter='latex'; c.FontSize=fs; 
  


%================================================================================       
% Text: Cumulative time-averaged flow. 
    axes('Position',[0.0  0.2  0.1 0.1]);     
    set(gcf,'color','w');axis equal; box on; box off; axis off;
    %text(0,0,0,{['Cumulative time-' char(10) 'averaged flow']},'FontName','times','FontSize',fs,'FontWeight','bold');   
    text(0,0,0,sprintf('Cumulative time- \n averaged flow'),'FontName','times','FontSize',fs,'FontWeight','bold');   
    
  
%% Save the video.
    drawnow;
    frame(i_nt)=getframe(gcf);
    writeVideo(v,frame(i_nt));       
end


