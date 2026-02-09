function rxm_ReconFlow_InstantFar_Plot

% This function plots the far-field instant flow around the reconstructed human sperm, relative to the comoving frame of reference.
% Figures:
% 1. 3D flow. 
% 2. 2D slicestream.
% 3. Sperm snapshot.
% 4. Colorbar.

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
    [~,np_beat_temp] = min(abs(t_nd_flow-2*pi));
    np_beat = np_beat_temp-1;
    % Select the instant moment.
    inst = 8.2; % According to the consecutive flagellar shapes within one beat cycle and the far-field flow features' switch (pusher/puller), 5 moments are selected: 4/5/6/7/8.2.
    nt0 = round(inst/5*np_beat);  % The plotted moment of the instant flow.

    xg=linspace(xl_field(1),xl_field(2),Nx_field);
    yg=linspace(yl_field(1),yl_field(2),Ny_field); 
    zg=linspace(zl_field(1),zl_field(2),Nz_field); 
    [Xg,Yg,Zg]=meshgrid(xg,yg,zg);  
    
    Ug_nt0 = reshape(Ug(:,nt0),size(Xg,1),size(Xg,2),size(Xg,3));
    Vg_nt0 = reshape(Vg(:,nt0),size(Xg,1),size(Xg,2),size(Xg,3));
    Wg_nt0 = reshape(Wg(:,nt0),size(Xg,1),size(Xg,2),size(Xg,3));
    vel_raw=(Ug_nt0.*Ug_nt0+Vg_nt0.*Vg_nt0+Wg_nt0.*Wg_nt0).^0.5;
    vel_nd=vel_raw.*2*pi; % Transform the velocity units into 'flagellum length/ beat cycle'.
    vel_dim = vel_nd.*scale_vel;  % in units of '\mu m/s'.
  
    
    %=====================================================================================================
    % Get sperm shape.
    load XNodes_sp8_NoBoundStokeslet.mat; 
    % For the datasets used for figure plot, rather than for generating video, the time points of flow data is coarser than those of sperm mobility.
    % Find suitable time moment for the corresponding instant flow.
    nt_XNodes=find(t_nd==t_nd_flow(end));
    nt0_XNodes=find(t_nd==t_nd_flow(nt0));
    Q=Nhh+Ns;   
    Xtail_raw = Xs(Nhh+1:Q,1:nt_XNodes);
    Ytail_raw = Xs(Q+Nhh+1:2*Q,1:nt_XNodes); 
    Ztail_raw = Xs(2*Q+Nhh+1:3*Q,1:nt_XNodes); 
    % Spatial smooth of the flagellar shape, to remove the noise arising from the raw experimental data being reconstructed.
    ppx = 0.999; % 'ppx=1' means interpolation without cubic spline.
    ppy = 0.999;
    ppz = 0.999;  
    for i_nt = 1:nt_XNodes
        s_temp = arclength(Xtail_raw(:,i_nt), Ytail_raw(:,i_nt), Ztail_raw(:,i_nt));
        Xtail_sm(:,i_nt) = fnval(csaps(s_temp,Xtail_raw(:,i_nt),ppx),s_temp);
        Ytail_sm(:,i_nt) = fnval(csaps(s_temp,Ytail_raw(:,i_nt),ppy),s_temp);
        Ztail_sm(:,i_nt) = fnval(csaps(s_temp,Ztail_raw(:,i_nt),ppz),s_temp);
    end   
    
    % Get the mean flagellar shape direction, based on which the instant sperm shape relative to the comoving frame can be obtained later by rotation.
    Xtail_0=Xtail_sm-repmat(Xtail_sm(1,:),Ns,1);
    Ytail_0=Ytail_sm-repmat(Ytail_sm(1,:),Ns,1);
    Ztail_0=Ztail_sm-repmat(Ztail_sm(1,:),Ns,1);   
    Xtail_mean= mean(Xtail_0,2); %Ns*1  %The mean flagellar shape relative to the lab frame.
    Ytail_mean= mean(Ytail_0,2);
    Ztail_mean= mean(Ztail_0,2);
    dir_new = [1, 0, 0];
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
    hold on; axis equal; axis off; 
    %//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    XL=[xl_field(1)-2 xl_field(2)]; YL=[yl_field(1) yl_field(2)]; ZL=[zl_field(1) zl_field(2)];  % for flow.
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
    
    
    if 0 % Plot 2D slice planes.
    % X planes.
    XL1_mat = 0; %[-4 -1.5 0 1 4];
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
    
    
    
    % Plot flow streamlines.
    % cross-sectional flow at y=0.
    [sx_3D,sy_3D,sz_3D] = meshgrid(xg(1:1:end),yg(21),zg(1:1:end));
    Generate_FlowStreamlines(Xg,Yg,Zg,Ug_nt0,Vg_nt0,Wg_nt0,sx_3D,sy_3D,sz_3D, vel_dim);
    % cross-sectional flow at z=0 (3/4 plane).
    [sx_3D,sy_3D,sz_3D] = meshgrid(xg(1:1:end),yg(22:1:end),zg(21));
    Generate_FlowStreamlines(Xg,Yg,Zg,Ug_nt0,Vg_nt0,Wg_nt0,sx_3D,sy_3D,sz_3D, vel_dim);
    [sx_3D,sy_3D,sz_3D] = meshgrid(xg(23:1:end),yg(1:1:20),zg(21));
    Generate_FlowStreamlines(Xg,Yg,Zg,Ug_nt0,Vg_nt0,Wg_nt0,sx_3D,sy_3D,sz_3D, vel_dim);
       
    shading interp
    colorbar; 
    caxis([0 0.1]); 
    
    
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if 1   % Arrows indicating the streamlines' directions, for full flow (reconstructed sperm, sp#8, no boundary).
        % arrows located at the top.
        [sx_cone1,sy_cone1,sz_cone1] = meshgrid(xg([23]),yg([3 26 40]),zg([4 18 36]));  
        sx_cone1([1 3],:,:)=NaN; sy_cone1([1 3],:,:)=NaN; sz_cone1([1 3],:,:)=NaN;
        sx_cone1(2,:,[1 2])=NaN; sy_cone1(2,:,[1 2])=NaN; sz_cone1(2,:,[1 2])=NaN;
        % arrows located at the bottom.
        [sx_cone2,sy_cone2,sz_cone2] = meshgrid(xg([23]),yg([3 26 40]),zg([6 18 41]));  
        sx_cone2([1 3],:,:)=NaN; sy_cone2([1 3],:,:)=NaN; sz_cone2([1 3],:,:)=NaN;
        sx_cone2(2,:,[2 3])=NaN; sy_cone2(2,:,[2 3])=NaN; sz_cone2(2,:,[2 3])=NaN;
        % arrows located at the left.
        [sx_cone3,sy_cone3,sz_cone3] = meshgrid(xg([10 23 40]),yg([2 26 36]),zg([22]));  
        sx_cone3(:,[1 3],:)=NaN; sy_cone3(:,[1 3],:)=NaN; sz_cone3(:,[1 3],:)=NaN;
        sx_cone3([1 2],2,:)=NaN; sy_cone3([1 2],2,:)=NaN; sz_cone3([1 2],2,:)=NaN;
       % arrows located at the right.
        %[sx_cone4,sy_cone4,sz_cone4] = meshgrid(xg([10 23 40]),yg([4 26 40]),zg([22]));  
        %sx_cone4(:,[1 3],:)=NaN; sy_cone4(:,[1 3],:)=NaN; sz_cone4(:,[1 3],:)=NaN;
        %sx_cone4([2 3],2,:)=NaN; sy_cone4([2 3],2,:)=NaN; sz_cone4([2 3],2,:)=NaN;
        % arrows located in the front.
        [sx_cone5,sy_cone5,sz_cone5] = meshgrid(xg([6 23 36]),yg([24]),zg([8 22 40]));  
        sx_cone5(:,:,[1 3])=NaN;  sy_cone5(:,:,[1 3])=NaN;   sz_cone5(:,:,[1 3])=NaN;
        sx_cone5(:,[1 2],2)=NaN;  sy_cone5(:,[1 2],2)=NaN;   sz_cone5(:,[1 2],2)=NaN;
        % arrows located at the back.
        [sx_cone6,sy_cone6,sz_cone6] = meshgrid(xg([2 23 36]),yg([24]),zg([8 22 40]));  
        sx_cone6(:,:,[1 3])=NaN;  sy_cone6(:,:,[1 3])=NaN;   sz_cone6(:,:,[1 3])=NaN;
        sx_cone6(:,[2 3],2)=NaN;  sy_cone6(:,[2 3],2)=NaN;   sz_cone6(:,[2 3],2)=NaN;
        
        Ug_cone = Ug_nt0./vel_raw;  Vg_cone = Vg_nt0./vel_raw;  Wg_cone = Wg_nt0./vel_raw;
        scal = 3; 
        Plot_StreamCone(fig,Xg,Yg,Zg,Ug_cone,Vg_cone,Wg_cone,sx_cone1,sy_cone1,sz_cone1,scal)
        Plot_StreamCone(fig,Xg,Yg,Zg,Ug_cone,Vg_cone,Wg_cone,sx_cone2,sy_cone2,sz_cone2,scal)
        Plot_StreamCone(fig,Xg,Yg,Zg,Ug_cone,Vg_cone,Wg_cone,sx_cone3,sy_cone3,sz_cone3,scal)
        %Plot_StreamCone(fig,Xg,Yg,Zg,Ug_cone,Vg_cone,Wg_cone,sx_cone4,sy_cone4,sz_cone4,scal)
        Plot_StreamCone(fig,Xg,Yg,Zg,Ug_cone,Vg_cone,Wg_cone,sx_cone5,sy_cone5,sz_cone5,scal)
        Plot_StreamCone(fig,Xg,Yg,Zg,Ug_cone,Vg_cone,Wg_cone,sx_cone6,sy_cone6,sz_cone6,scal)
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




    



%% Fig. 3: Sperm snapshot.

if 1

% Select a flow-plot moment.
    load FarFlow_CF_sp8_NoBound.mat;
    [~,np_beat_temp] = min(abs(t_nd_flow-2*pi));
    np_beat = np_beat_temp-1;
    % Select the instant moment.
    inst = 8.2; %4; 5; 6; 7; 8.2;
    %nbc = length(t_nd_flow)/np_beat;
    nt0 = round(inst/5*np_beat);  % The plotted moment of the instant flow.
    % Transform the non-dimensional time into dimensional form.
    sp=8;
    load FreeSperm_Frequency.mat; % read the beat frequency (HF WF frequency) to encode the WF color
    Freq_WF = HF_Freq(sp);
    t_dim_flow = t_nd_flow/(2*pi)/Freq_WF;
    clearvars -except sp t_inst_dim t_dim_flow nt0
    
    
% Process the waveform data.
load 03_Head_Fixed_Frame_Method_Rotations_Traces_raw_5_Microns_Difraction_head_spin_rxmSSppxyz.mat;
Xtail = X{sp}; 
Ytail = Y{sp}; 
Ztail = Z{sp}; 
nt=size(Xtail,2);
ns=size(Xtail,1);

for i_nt=1:nt
    arc_mat(:,i_nt) = arclength(Xtail(:,i_nt), Ytail(:,i_nt), Ztail(:,i_nt));
end
arc_mean = mean(arc_mat(end,:)); % measured as 46.35.

% Temporally interpolate the original waveform.
dt=1/90; %sampling frequency of the experimental imaging.
tRange=0:dt:(nt-1)*dt;
t0 = linspace(0,tRange(end),nt);
t1 = t_dim_flow; %t1=linspace(0,tRange(end),5*nt);
s=1:1:ns;
Xtail = interp2(t0,s',Xtail,t1,s'); 
Ytail = interp2(t0,s',Ytail,t1,s'); 
Ztail = interp2(t0,s',Ztail,t1,s'); 
nt=size(Xtail,2);
%t_inst_dim = t_dim_flow(nt0);

% Generate sperm head: ideal sperm model.
[Xhead_up,Yhead_up,Zhead_up, Xhead_down,Yhead_down,Zhead_down] = Get_BF_IdealSperm_Head(arc_mean);  %Xhead_up: M1*M2.


fs = 36;
figure(5)
set(gcf,'color','w','units','normalized','outerposition',[0 0 1 1]) 
clf; hold on;
box on; axis equal; axis off;
view(-45,16); 
XL = [-6 max(Xtail(:))];  YL = [min(Ytail(:)) max(Ytail(:))];  ZL = [min(Ztail(:)) max(Ztail(:))]; 
xlim([XL(1),XL(2)]);ylim([YL(1),YL(2)]);zlim([ZL(1),ZL(2)]);
xlh=xlabel('$\xi_3 (\mu m)$ ','interpreter','latex','FontSize',fs);
ylh=ylabel('$\xi_1 (\mu m)$ ','interpreter','latex','FontSize',fs);
zlh=zlabel('$\xi_2 (\mu m)$ ','interpreter','latex','FontSize',fs,'Rotation',0);
% Plot the head.
source_ligth = [90 80];
k = [0.9 1 1 5]; 
s1= surfl(Xhead_up,Yhead_up,Zhead_up,source_ligth,k );
s2= surfl(Xhead_down,Yhead_down,Zhead_down,source_ligth,k );
s1.EdgeColor='none' ;
s1.FaceColor=  [0.6 0 0 ];
s2.EdgeColor='none' ;
s2.FaceColor=  [0 0.6 0]; 
camlight right
lighting gouraud  
% Plot tail waveform
for i_nt = 1:1:nt
   p = plot3(Xtail(:,i_nt),Ytail(:,i_nt),Ztail(:,i_nt),'color',[36 100 171]./255,'linewidth',2);
   % transparency
   p.Color(4) = 0.1;
end
plot3(Xtail(:,nt0),Ytail(:,nt0),Ztail(:,nt0),'color',[36 100 171]./255,'linewidth',5)

return
end

%% Fig. 4: Colorbar.

figure(4)
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


