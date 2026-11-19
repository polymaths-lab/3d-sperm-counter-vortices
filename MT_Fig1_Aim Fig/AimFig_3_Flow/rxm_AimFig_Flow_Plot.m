function rxm_AimFig_Flow_Plot

% This function plots the predicted 3D avergae flow around the observed human sperm, figure used for the 'aim figure' in the manuscript. 
% The flow field is relative to the comoving frame, and is averaged within complete trajectory revolution periods. 
% The aim figure is for visualization purpose only, so the sperm body in this function is not to exact scale.

clc
%sp=8;    

    
%% Step 1: determine the sperm body position.
    
    %============================================================================================================================
    % generate sperm head: humanoid size
    % creating mesh for sperm head provide by Hermes
    Ndiv = 60;
    radius_plot = 1.2/56;%normalized
    head_type = 1;    
    [Xhead,Yhead,Zhead] = sperm_head3D_shape(radius_plot, head_type, Ndiv);%size(Zhead)=[61 61]
    scal=1.4;  %1.7;%scal=80;%scal=56;   %This 'scal' parameter can be adjusted according to visualization requirement.
    Xhead=Xhead*scal; Yhead=Yhead*scal; Zhead=Zhead*scal;
    M1=size(Xhead,1);M2=size(Xhead,2);M=M1*M2;
    [~,Ind_Neck] = max(Xhead(:));
    % head upside
    [zup_row,zup_col]=find(Zhead>=0); %length(zup_row)=1891
    Zhead_up=nan(size(Zhead));Xhead_up=nan(size(Xhead));Yhead_up=nan(size(Yhead));
    for i = 1:length(zup_row)
        Zhead_up(zup_row(i),zup_col(i))=Zhead(zup_row(i),zup_col(i));
        Xhead_up(zup_row(i),zup_col(i))=Xhead(zup_row(i),zup_col(i));
        Yhead_up(zup_row(i),zup_col(i))=Yhead(zup_row(i),zup_col(i));
    end
    % head downside
    [zdown_row,zdown_col]=find(Zhead<=0);%length(zdown_row)=1891
    Zhead_down=nan(size(Zhead));Xhead_down=nan(size(Xhead));Yhead_down=nan(size(Yhead));
    for i = 1:length(zdown_row)
        Zhead_down(zdown_row(i),zdown_col(i))=Zhead(zdown_row(i),zdown_col(i));
        Xhead_down(zdown_row(i),zdown_col(i))=Xhead(zdown_row(i),zdown_col(i));
        Yhead_down(zdown_row(i),zdown_col(i))=Yhead(zdown_row(i),zdown_col(i));
    end           
    % Rotate the head to align with the average flow field.
    b1T = [-1; 0; 0];
    b2N = [0; -1; 0];
    b3B = [0; 0; 1];
    B = [b1T b2N b3B];
    Head_up = ApplyRotationMatrix(B,[Xhead_up(:);Yhead_up(:);Zhead_up(:)]);
    [Xhead_up, Yhead_up, Zhead_up] = ExtractComponents(Head_up);
    Xhead_up = reshape(Xhead_up,M1,M2);
    Yhead_up = reshape(Yhead_up,M1,M2);
    Zhead_up = reshape(Zhead_up,M1,M2);
    Head_down = ApplyRotationMatrix(B,[Xhead_down(:);Yhead_down(:);Zhead_down(:)]);
    [Xhead_down, Yhead_down, Zhead_down] =ExtractComponents(Head_down);
    Xhead_down = reshape(Xhead_down,M1,M2);
    Yhead_down = reshape(Yhead_down,M1,M2);
    Zhead_down = reshape(Zhead_down,M1,M2);
    
        
    
    %============================================================================================================================
    % get flagellum shape.
    load XNodes_sp8_NoBoundStokeslet.mat;
    load NearFlow_CF_sp8_NoBound.mat;
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
    Xtail_1=Xtail_0-repmat(Xtail_0(1,:),Ns,1);
    Ytail_1=Ytail_0-repmat(Ytail_0(1,:),Ns,1);
    Ztail_1=Ztail_0-repmat(Ztail_0(1,:),Ns,1);    
    xtail_mean_temp= mean(Xtail_1,2); %Ns*1 
    ytail_mean_temp= mean(Ytail_1,2);
    ztail_mean_temp= mean(Ztail_1,2);
    % Align flagellum to get the comoving waveform, that's also how we get
    % the aligned flow field box, so that the average flagellum shape can
    % be generally parallel to the X axis of the aligned field box.
    dir_new=[1 0 0];
    dir_raw=[xtail_mean_temp(1)-xtail_mean_temp(end),ytail_mean_temp(1)-ytail_mean_temp(end),ztail_mean_temp(1)-ztail_mean_temp(end)];
    nt_plot = 33; % choose a plot moment for the flagellum.
    [Xtail_plot,Ytail_plot,Ztail_plot] = rxm_dir_align(dir_raw,dir_new,Xtail_1(:,nt_plot),Ytail_1(:,nt_plot),Ztail_1(:,nt_plot));
    Xtail_plot=Xtail_plot+Xhead(Ind_Neck);
    Ytail_plot=Ytail_plot+Yhead(Ind_Neck);
    Ztail_plot=Ztail_plot+Zhead(Ind_Neck);
    
  

    
%% Step 2: get the average flow velocity.

    % 3D grid space of the local flow field.
    xg=linspace(xl_field(1),xl_field(2),Nx_field);
    yg=linspace(yl_field(1),yl_field(2),Ny_field); 
    zg=linspace(zl_field(1),zl_field(2),Nz_field); 
    [Xg,Yg,Zg]=meshgrid(xg,yg,zg);   
    
    % Calculate the average flow velocity.
    Ug_mean = mean(Ug,2);  
    Vg_mean = mean(Vg,2);  
    Wg_mean = mean(Wg,2); 
    Ug_mean = reshape(Ug_mean,size(Xg,1),size(Xg,2),size(Xg,3));
    Vg_mean = reshape(Vg_mean,size(Xg,1),size(Xg,2),size(Xg,3));
    Wg_mean = reshape(Wg_mean,size(Xg,1),size(Xg,2),size(Xg,3));
    % Velocity units: flagellum length/ beat cycle (non-dimensional).
    vel=(Ug_mean.*Ug_mean+Vg_mean.*Vg_mean+Wg_mean.*Wg_mean).^0.5; 
    vel=vel.*2*pi;  

    
    
    
    
%% Step 3: plot the time-averaged 3D streamlines.

fig = 3;
figure(fig)
clf; hold on; axis equal; 
% For exp. sperm.
view([120 10]);  
set(gcf,'color','w');
%ax = gca; ax.FontSize =16;  ax.TickLabelInterpreter = 'latex';
%ax = gca; ax.FontSize = 30;  ax.TickLabelInterpreter = 'latex';
%xlh=xlabel('$x$','interpreter','latex','FontSize',40);
%ylh=ylabel('$y$','interpreter','latex','FontSize',40);
%zlh=zlabel('$z$','interpreter','latex','FontSize',40,'Rotation',0);
axis off; 


%============================================================================================================================
% Plot sperm body.
PlotTailCylinder_AimFig(Xtail_plot,Ytail_plot,Ztail_plot,1,fig)
CO_up(:,:,1) = 177/255*ones(size(Xhead_up)); CO_up(:,:,2) = 24/255*ones(size(Xhead_up)); CO_up(:,:,3) = 45/255*ones(size(Xhead_up));
CO_down(:,:,1) = 218/255*ones(size(Xhead_down)); CO_down(:,:,2) = 207/255*ones(size(Xhead_down)); CO_down(:,:,3) = 168/255*ones(size(Xhead_down));
s1= surf(Xhead_up,Yhead_up,Zhead_up,CO_up);
s2= surf(Xhead_down,Yhead_down,Zhead_down,CO_down);
s1.EdgeColor='none' ;
s2.EdgeColor='none' ;
                
            
             
%============================================================================================================================
% Plot the 3D streamlines. For better visualization, we need to assign more than one groups of streamlines.  

% Streamline Group 1: near the flagellum tip.
for i_g1 = 1:4
    
    % 1.Indicate the start of the 3D streamlines.
    if i_g1==1
        [sx_3D,sy_3D,sz_3D] = meshgrid(xg(10:11),yg(11),zg(10)); 
    elseif i_g1==2
        [sx_3D,sy_3D,sz_3D] = meshgrid(xg(12),yg(10:11),zg(8:2:10)); 
    elseif i_g1==3
        [sx_3D,sy_3D,sz_3D] = meshgrid(xg(13),yg(12:13),zg(8));  
    elseif i_g1==4
        [sx_3D,sy_3D,sz_3D] = meshgrid(xg(14),yg(14),zg(8));  
    end
    
    % 2. Locate each vertices along the 3D streamlines.        
    vert = stream3(Xg,Yg,Zg,Ug_mean,Vg_mean,Wg_mean,sx_3D,sy_3D,sz_3D);
    % Get the velocity values along streamlines.
    % 'Vetices' is a cell array and each cell in this case has size '?*3', which gives the [X Y Z] coordinates along each streamline.
    % That's why we need to interpolate the velocity field to get the velocity at these vertices.
    vert_mat = [0 0 0];
    nsl = length(vert);  %number of streamlines
    nvert_sl = nan(nsl,1); % record the number of vertices for each streamline.
    for i_sl = 1:nsl
        vert_temp = vert{i_sl};
        nvert_sl(i_sl) = size(vert_temp,1);
        vert_mat = [vert_mat; vert_temp];
    end
    vert_mat = vert_mat(2:end,:);
    
    % 3. Map streamline colors according to velocity values.
    Vq = interp3(Xg,Yg,Zg,vel,vert_mat(:,1),vert_mat(:,2),vert_mat(:,3));
    nv = 0;
    for i_sl = 1:nsl
        nv0 = nv+1;
        nv = nv+nvert_sl(i_sl);
       vert_temp = { vert{i_sl} };
        h = streamtube(vert_temp ,0.01); 
        h.CData = repmat(Vq(nv0:nv),1,21);   
    end
    shading interp
    colormap(gca,'parula'); caxis([0 0.01]);
end




% Streamline Group 2: near sperm head.
for i_g2=1:2
    
    % 1.Indicate the start of the 3D streamlines.   
    if i_g2==1
        [sx_3D,sy_3D,sz_3D] = meshgrid(xg(13),yg(9:10),zg(8:9)); 
    elseif i_g2==2
        [sx_3D,sy_3D,sz_3D] = meshgrid(xg(14:15),yg(13),zg(8:10)); 
    end
    
    % 2. Locate each vertices along the 3D streamlines.
    vert = stream3(Xg,Yg,Zg,Ug_mean,Vg_mean,Wg_mean,sx_3D,sy_3D,sz_3D);
    vert_mat = [0 0 0];
    nsl = length(vert);  %number of streamlines
    nvert_sl = nan(nsl,1); % record the number of vertices for each streamline.
    for i_sl = 1:nsl
        vert_temp = vert{i_sl};
        nvert_sl(i_sl) = size(vert_temp,1);
        vert_mat = [vert_mat; vert_temp];
    end
    vert_mat = vert_mat(2:end,:);

    % 3. Map streamline colors according to velocity values.
    Vq = interp3(Xg,Yg,Zg,vel,vert_mat(:,1),vert_mat(:,2),vert_mat(:,3));
    nv = 0;
    for i_sl = 1:nsl
        nv0 = nv+1;
        nv = nv+nvert_sl(i_sl);
        vert_temp = { vert{i_sl} };
        h = streamtube(vert_temp ,0.01); 
        h.CData = repmat(Vq(nv0:nv),1,21);   
    end
    shading interp
    colormap(gca,'parula'); caxis([0 0.01]);
end





XL=[-1.2 0.4]; YL=[-0.22 0.3]; ZL=[-0.22 0.2];  %NB, BH0.2
xlim([XL(1) XL(2)]); ylim([YL(1) YL(2)]); zlim([ZL(1) ZL(2)]); 

%camlight right
%lighting gouraud
l1=light;
l1.Color = [1 1 1];
l1.Style = 'local';  %'infinite';  %
l1.Position = [0.2 0.3 0.1];
l2=light;
l2.Color = [1 1 1];
l2.Style = 'local';  %'infinite';  %
l2.Position = [-0.8 0.3 0];
l3=light;
l3.Color = [1 1 1];
l3.Style = 'local';  %'infinite';  %
l3.Position = [0.3 -0.3 0];
lighting gouraud   %phong  %flat  %
material metal  %shiny     
