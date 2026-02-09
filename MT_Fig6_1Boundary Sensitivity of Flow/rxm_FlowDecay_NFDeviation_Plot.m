function rxm_FlowDecay_NFDeviation_Plot

clc


sp=8;

load 01_lab_frame_raw_traces_2017_2018_rxmSSppxyz.mat; 
arc_mat = arclength(X{sp},Y{sp},Z{sp});
arc_mean = mean(arc_mat(end,:));

load FreeSperm_Frequency.mat;  
scale_vel = arc_mean*HF_Freq(sp); % in units o '\mu m/s'.
 
clearvars -except scale_vel arc_mean

%% Get flow field.======================================================
 
load NearFlow_CF_sp8_NoBound.mat; 
% plotting grid to form the local flow field
    xg=linspace(xl_field(1),xl_field(2),Nx_field);
    yg=linspace(yl_field(1),yl_field(2),Ny_field); 
    zg=linspace(zl_field(1),zl_field(2),Nz_field); 
    [Xg,Yg,Zg]=meshgrid(xg,yg,zg); 
    Xg=Xg.*arc_mean; Yg=Yg.*arc_mean;  Zg=Zg.*arc_mean;  

% Calculate the average flow velocity of the moving
% field.===============================================================
    nt_flow=length(t_nd_flow);
    Ug_mean = mean(Ug(:,1:nt_flow),2);
    Vg_mean = mean(Vg(:,1:nt_flow),2);
    Wg_mean = mean(Wg(:,1:nt_flow),2);
    Ug_mean = reshape(Ug_mean,size(Xg,1),size(Xg,2),size(Xg,3));
    Vg_mean = reshape(Vg_mean,size(Xg,1),size(Xg,2),size(Xg,3));
    Wg_mean = reshape(Wg_mean,size(Xg,1),size(Xg,2),size(Xg,3));
    % Transform velocity units into 'flagellum length/ beat cycle',
    vel_NB=(Ug_mean.*Ug_mean+Vg_mean.*Vg_mean+Wg_mean.*Wg_mean).^0.5; 
    vel_NB=vel_NB.*2*pi;  %in units of 'flagellum length/ beat cycle',
    vel_NB = vel_NB.*scale_vel;  % in units of '\mu m/s'.
    clear Ug Vg Wg

    
load NearFlow_CF_sp8_BoundH0.2.mat; 
% Calculate the average flow velocity of the moving
% field.===============================================================
    nt_flow=length(t_nd_flow);
    Ug_mean = mean(Ug(:,1:nt_flow),2);
    Vg_mean = mean(Vg(:,1:nt_flow),2);
    Wg_mean = mean(Wg(:,1:nt_flow),2);
    Ug_mean = reshape(Ug_mean,size(Xg,1),size(Xg,2),size(Xg,3));
    Vg_mean = reshape(Vg_mean,size(Xg,1),size(Xg,2),size(Xg,3));
    Wg_mean = reshape(Wg_mean,size(Xg,1),size(Xg,2),size(Xg,3));
    % Transform velocity units into 'flagellum length/ beat cycle',
    vel_BH=(Ug_mean.*Ug_mean+Vg_mean.*Vg_mean+Wg_mean.*Wg_mean).^0.5; 
    vel_BH=vel_BH.*2*pi;  %in units of 'flagellum length/ beat cycle',
    vel_BH = vel_BH.*scale_vel;  % in units of '\mu m/s'.
    
    
    vel_dev = (vel_BH-vel_NB)./vel_NB;  


%% Plot the velocity (deviation) decay.    
    
    if 0
    [Mvel,Ind_Mvel] = max(vel(:));   
    [Iy,Ix,Iz] = ind2sub(size(Ug_mean),Ind_Mvel);
    [Mvel,Xg(Iy,Ix,Iz),Yg(Iy,Ix,Iz),Zg(Iy,Ix,Iz)]
    end
    % Locate the center in space for decay measure.
    % The far/near-field distance is measured from the middle of flagellum,
    % so that we have a 'near volume' with the radius of half flagellum
    % length that can basically cover the entire sperm body.
    [MinX,Ix] = min(abs((xg+0.45))); % for exp. sperm.
    %[MinX,Ix] = min(abs((xg-0.45))); % for num. sperm.
    [MinY,Iy] = min(abs((yg-0)));
    [MinZ,Iz] = min(abs((zg-0)));
    [MinX*arc_mean, MinY*arc_mean, MinZ*arc_mean, Xg(Iy,Ix,Iz),Yg(Iy,Ix,Iz),Zg(Iy,Ix,Iz)]
    
    
    fs=120;
    figure(2)
    clf; hold on; set(gcf,'color','w'); box on;
    daspect([1 0.06 1]); %axis square;
    ax = gca;
    ax.FontSize = fs; 
    ax.TickLabelInterpreter = 'latex';
    xlabel('$r \ (\mu m)$','interpreter','latex','FontSize',fs+20);
    %ylabel('$\Delta u \ (\mu m/s)$','interpreter','latex','FontSize',fs+10); % absolute deviation
    ylabel('$\Delta u \ (\%)$','interpreter','latex','FontSize',fs+20); % relative deviation
    %xlim([25 550]*60/50); ylim([0.000001 50]); % for far-field, exp. sperm (sp#8).
    %xlim([5 60]); ylim([0.01 30]); % for near-field. exp. sperm (sp#8, 23)
    %xlim([0 60]); ylim([-1 1.7]); % for decay deviation of near-field, in units of um/s. exp. sperm (sp#8)
    xlim([0 46]); ylim([-0.75 0.75]); xticks([0 20 40]); yticks([-0.5:0.5:0.5]); yticklabels({'-50','0','50'});% for decay deviation of near-field, in units of %. exp. sperm (sp#8)
    set(gca,'linewidth',1);

    % Velocity decay along axis X.
    for i_Indx=1:Ix
        xg1(i_Indx)=Xg(Iy,Ix,Iz)-Xg(Iy,Ix-i_Indx+1,Iz); 
        vel_xg1(i_Indx)=vel_dev(Iy,Ix-i_Indx+1,Iz);
    end        
    plot(xg1(:),vel_xg1(:),':','color',[139 0 139]/255,'linewidth',10);  %purple
    xg2=Xg(Iy,Ix:end,Iz)-Xg(Iy,Ix,Iz); vel_xg2=vel_dev(Iy,Ix:end,Iz);
    plot(xg2(:),vel_xg2(:),'color',[139 0 139]/255,'linewidth',10);  %purple
    
    % Velocity decay along axis Y.
    for i_Indy=1:Iy
        yg1(i_Indy)=Yg(Iy,Ix,Iz)-Yg(Iy-i_Indy+1,Ix,Iz); 
        vel_yg1(i_Indy)=vel_dev(Iy-i_Indy+1,Ix,Iz);
    end
    plot(yg1(:),vel_yg1(:),':','color',[255 215 0]/255,'linewidth',10);  %yellow/ golden  
    yg2=Yg(Iy:end,Ix,Iz); vel_yg2=vel_dev(Iy:end,Ix,Iz);
    plot(yg2(:),vel_yg2(:),'color',[255 215 0]/255,'linewidth',10);  %yellow/ golden  
        
    % Velocity decay along axis Z.
    for i_Indz=1:Iz 
        zg1(i_Indz)=Zg(Iy,Ix,Iz)-Zg(Iy,Ix,Iz-i_Indz+1); 
        vel_zg1(i_Indz)=vel_dev(Iy,Ix,Iz-i_Indz+1);
    end
    % for cases with 1 boundary underneath, the -Z direction points
    % seleted to be plotted vary with different boundary-height cases.
    % for near-field:
     plot(zg1(1:3),vel_zg1(1:3),':', 'color',[0 128 0]/255,'linewidth',10);  %green
    % plot(zg1(1:6),vel_zg1(1:6),'c:','linewidth',8)  % boundary height=0.5
    %plot(zg1(:),vel_zg1(:),'c:','linewidth',8) % boundary height=1, and for cases without a boundary.
    % for far-field: we don't plot -Z direction when with boundary, because the sparse spatial
    % discretisation barely records the close distance to the boundary.
    zg2=Zg(Iy,Ix,Iz:end); vel_zg2=vel_dev(Iy,Ix,Iz:end);
    plot(zg2(:),vel_zg2(:), 'color',[0 128 0]/255,'linewidth',10);  %green
       

    if 0
        % Standard line, with r^{-1}/r^{-2}/r^{-3} velocity decay.
    if 1  % for far-field.
        plot([25 2500],[3 0.0003],'k','linewidth',4);    % For numerical model alpha=0.6 and experimental sperm8.
        text(100,1,'$r^{-2}$','interpreter','latex','fontsize',fs)
        plot([25 2500],[0.1 0.0000001],'k','linewidth',4);    % For numerical model alpha=0.6 and experimental sperm8.
        text(35,0.001,'$r^{-3}$','interpreter','latex','fontsize',fs)
    else  % for near-field
        plot([5 50],[30 3],'k','linewidth',4);    % For experimental sperm8.
        text(14,15,'$r^{-1}$','interpreter','latex','fontsize',fs)
        plot([5 50],[2 0.02],'k','linewidth',4);    % For experimental sperm8.
        text(7,0.2,'$r^{-2}$','interpreter','latex','fontsize',fs)   
    end
    end
    
    % Standard line, with r^{-1} velocity decay.
    %plot([10 1000],[10 0.1],'k:','linewidth',2);    % For numerical model, alpha=0.4
    %plot([10 1000],[4 0.04],'k:','linewidth',2);    % For experimental sperm8.
    %text(10,4,'$r^{-1}$','interpreter','latex','fontsize',30)
    
    %text(35,0.00001,'$H=\infty$','interpreter','latex','fontsize',fs-10)
    %text(35,0.00001,'$H=1$','interpreter','latex','fontsize',fs-10)
    % text position for near-field plot: [6,0.02]
    % text position for far-field plot: [35,0.00001]

