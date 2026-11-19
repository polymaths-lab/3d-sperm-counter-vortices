function rxm_IdealFlow_InstantFar_PushPullSwitch_Plot

% This function plots the instant far-field flows of the idealized
% model: pusher-I, II, III and puller-I, II, III.

clc


%% Plot: XYZ coordinate system/ colorbar.
%================================================================================
%================================================================================
if 0
    fig = 1; figure(fig) 
    clf;    
    hold on; set(gcf,'color','w'); axis equal; box on; box off; axis off;
    view([-50 15]); 
    xlabel('x'); ylabel('y'); zlabel('z');
if 0    
    % Plot Cartesian axes.
    fs = 32;  
    O =[1 0 0];  Cx=[-1; 0; 0];  Cy=[0; -1; 0];  Cz=[0; 0; 1];  C=[Cx Cy Cz];
    Plot_Cartesian_Axes_V1(fig,O,C) 
    camlight right
    lighting gouraud
    % Text: axes labels.
    %text(0.4,0,0,{'$x_c$'},'FontName','times','interpreter','latex','FontSize',fs,'FontWeight','bold');   
    %text(1,-0.43,0,{'$y_c$'},'FontName','times','interpreter','latex','FontSize',fs,'FontWeight','bold');   
    %text(1,0.1,0.52,{'$z_c$'},'FontName','times','interpreter','latex','FontSize',fs,'FontWeight','bold'); 
elseif 1
    % Plot the colorbar.
    fs = 42;  
    caxis([0 0.0001]);      
    c=colorbar; 
    c.Position=[0.67 0.35 0.008 0.57];  
    c.Label.String = 'Flow velocity'; 
    %text(15,-17.5,-4,{'Flow velocity'},'FontName','times','FontSize',fs,'rotation',90);  %'FontWeight','bold'  
    c.Ticks=[0 0.0001]; c.TickLabels={'0','$10^{-4}$'}; 
    c.FontName = 'times'; c.TickLabelInterpreter='latex'; c.FontSize=fs; 
end

return
end
    

    






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

t_nd_flow = tRange_flow/(2*pi);
instant = roundn(t_nd_flow,-3);
i_nt = find(instant==0.78);  %instant=[0.56 0.15 0.48 0.33 0.4 0.78]
    

% Main figure: 3D comoving-frame flow.
%================================================================================
%================================================================================
    %Get instant flow.
    Ug_nt0 = reshape(Ug(:,i_nt),size(Xg,1),size(Xg,2),size(Xg,3));
    Vg_nt0 = reshape(Vg(:,i_nt),size(Xg,1),size(Xg,2),size(Xg,3));
    Wg_nt0 = reshape(Wg(:,i_nt),size(Xg,1),size(Xg,2),size(Xg,3));
    vel_raw=(Ug_nt0.*Ug_nt0+Vg_nt0.*Vg_nt0+Wg_nt0.*Wg_nt0).^0.5;
    vel_nd=vel_raw.*2*pi; % Transform the velocity units into 'flagellum length/ beat cycle'.    
       
    % Plot the 3D comoving-frame flow.
    clf;    
    %set(gca,'position',[0.09 0.3 0.7 0.7]);  %get(gca,'position')    
    hold on; set(gcf,'color','w'); axis equal; box on; box off; axis off;
    view([-50 15]); %view([-35 20]);  % for idealized model's far average flow. 
    %XL=[xl_field(1)-2.2 xl_field(2)+2.6]; YL=[yl_field(1)-1.7 yl_field(2)+1.8]; ZL=[zl_field(1)-2.2 zl_field(2)+2.2];
    XL=[xl_field(1)-14 xl_field(2)+14]; YL=[yl_field(1)-13 yl_field(2)+14]; ZL=[zl_field(1)-7 zl_field(2)+8];
    xlim([XL(1) XL(2)]); ylim([YL(1) YL(2)]); zlim([ZL(1) ZL(2)]);  
    xlabel('x'); ylabel('y'); zlabel('z');
  
    % Plot flow streamlines.
    % cross-sectional flow at y=0.
    [sx_3D,sy_3D,sz_3D] = meshgrid(xg(1:2:end),yg(21),zg(1:2:end)); 
    %[sx_3D,sy_3D,sz_3D] = meshgrid(xg(1:10:end),yg(21),zg(1:10:end)); 
    Generate_FlowStreamlines(Xg,Yg,Zg,Ug_nt0,Vg_nt0,Wg_nt0,sx_3D,sy_3D,sz_3D, vel_nd);       
    % cross-sectional flow at z=0.
    [sx_3D,sy_3D,sz_3D] = meshgrid(xg(1:2:end),yg(1:2:end),zg(21)); 
    %[sx_3D,sy_3D,sz_3D] = meshgrid(xg(1:10:end),yg(1:10:end),zg(21)); 
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
    Plot_StreamCone(fig,Xg,Yg,Zg+3,Ug_cone,Vg_cone,Wg_cone,sx_cone1,sy_cone1,sz_cone1+3,1.7*scal) % arrows located at the top.    
    Plot_StreamCone(fig,Xg,Yg,Zg-3,Ug_cone,Vg_cone,Wg_cone,sx_cone2,sy_cone2,sz_cone2-3,1.7*scal) % arrows located at the bottom.
    Plot_StreamCone(fig,Xg,Yg-8,Zg,Ug_cone,Vg_cone,Wg_cone,sx_cone3,sy_cone3-8,sz_cone3,2*scal) % arrows located at the left.
    Plot_StreamCone(fig,Xg,Yg+9,Zg,Ug_cone,Vg_cone,Wg_cone,sx_cone4,sy_cone4+9,sz_cone4,2*scal) % arrows located at the right.
    Plot_StreamCone(fig,Xg-8,Yg,Zg,Ug_cone,Vg_cone,Wg_cone,sx_cone5-8,sy_cone5,sz_cone5,2*scal) % arrows located in the front.
    Plot_StreamCone(fig,Xg+8,Yg,Zg,Ug_cone,Vg_cone,Wg_cone,sx_cone6+8,sy_cone6,sz_cone6,2*scal) % arrows located at the back.
    camlight right
    lighting gouraud
    
    caxis([0 0.0001]);      
    if 0
    % Plot the colorbar.
    c=colorbar; 
    c.Position=[0.67 0.35 0.005 0.57];  
    %c.Label.String = 'Flow velocity'; 
    text(15,-17.5,-4,{'Flow velocity'},'FontName','times','FontSize',fs,'rotation',90);  %'FontWeight','bold'  
    c.Ticks=[0 0.0001]; c.TickLabels={'0','$10^{-4}$'}; 
    c.FontName = 'times'; c.TickLabelInterpreter='latex'; c.FontSize=fs; 
    end
    
    % Plot 2D slice planes.
    if 1 
        if 0
        % X planes.
        XL1_mat = 0.5;  %[...]
        for i_x = 1:length(XL1_mat)
            XL_temp=[XL1_mat(i_x) XL(2)];
            YL_temp=[YL(1) YL(2)];
            ZL_temp=[ZL(1) ZL(2)]; 
            [xgrid,ygrid,zgrid] = meshgrid(linspace(XL_temp(1),XL_temp(2),10),...
                    linspace(YL_temp(1),YL_temp(2),10),...
                    linspace(ZL_temp(1),ZL_temp(2),10));     
            %CO(:,:,1) = ones(10)*.0; % red 
            %CO(:,:,2) = ones(10)*.4470; % green 
            %CO(:,:,3) = ones(10)*0.7410; % blue 
            CO(:,:,1) = ones(10)*1; % red 
            CO(:,:,2) = ones(10)*0; % green 
            CO(:,:,3) = ones(10)*0; % blue 
            surfX = surf(reshape(xgrid(:,1,:),[10 10]),reshape(ygrid(:,1,:),[10 10]),reshape(zgrid(:,1,:),[10,10]), CO,'FaceAlpha',0.1);  
            surfX.EdgeColor='none';
        end
        end
        if 1
        % Y plane.
        XL_temp=[xl_field(1)-10 xl_field(2)+10];  %XL_temp=[XL(1) XL(2)];
        YL_temp=[0 yl_field(2)]; %YL_temp=[0 YL(2)];
        ZL_temp=[zl_field(1)-5 zl_field(2)+5];   %ZL_temp=[ZL(1) ZL(2)]; 
        [xgrid,ygrid,zgrid] = meshgrid(linspace(XL_temp(1),XL_temp(2),10),...
                    linspace(YL_temp(1),YL_temp(2),10),...
                    linspace(ZL_temp(1),ZL_temp(2),10));     
        %CO(:,:,1) = ones(10)*.0; % red 
        %CO(:,:,2) = ones(10)*.4470; % green 
        %CO(:,:,3) = ones(10)*0.7410; % blue 
        CO(:,:,1) = ones(10)*1; % red 
        CO(:,:,2) = ones(10)*0; % green 
        CO(:,:,3) = ones(10)*0; % blue 
        surfY = surf(reshape(xgrid(1,:,:),[10 10]),reshape(ygrid(1,:,:),[10 10]),reshape(zgrid(1,:,:),[10,10]), CO,'FaceAlpha',0.1);  
        surfY.EdgeColor='none';
        end
        if 0
        % Z plane.          
        XL_temp=[xl_field(1)-10 xl_field(2)+10];  %XL_temp=[XL(1) XL(2)];
        YL_temp=[yl_field(1)-10 yl_field(2)+10];  %YL_temp=[YL(1) YL(2)];
        ZL_temp=[0 zl_field(2)];  %ZL_temp=[0 ZL(2)]; 
        [xgrid,ygrid,zgrid] = meshgrid(linspace(XL_temp(1),XL_temp(2),10),...
                    linspace(YL_temp(1),YL_temp(2),10),...
                    linspace(ZL_temp(1),ZL_temp(2),10));     
        %CO(:,:,1) = ones(10)*.0; % red 
        %CO(:,:,2) = ones(10)*.4470; % green 
        %CO(:,:,3) = ones(10)*0.7410; % blue 
        CO(:,:,1) = ones(10)*1; % red 
        CO(:,:,2) = ones(10)*0; % green 
        CO(:,:,3) = ones(10)*0; % blue 
        surfZ = surf(reshape(xgrid(:,:,1),[10 10]),reshape(ygrid(:,:,1),[10 10]),reshape(zgrid(:,:,1),[10,10]), CO,'FaceAlpha',0.1);  
        surfZ.EdgeColor='none';
        end
    end


