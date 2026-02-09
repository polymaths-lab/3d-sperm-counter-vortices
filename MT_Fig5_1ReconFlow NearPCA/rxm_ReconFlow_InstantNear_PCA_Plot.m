function rxm_ReconFlow_InstantNear_PCA_Plot

% This function plots the instant flow field, as well as its PCA reconstruction, around the reconstructed human sperm, relative to the comoving frame of reference.
% Figures:
% 1. 3D flow (raw and PCA-recon., and PCA modes).
% 2. PCA cumulative variance.
% 3. PCA time coefficients.
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

if 0
    op = 1;
    % op=1: plot the raw instant flow.
    % op=2: plot the PCA-reconstructed instant flow.
    % op=3: plot the PCA mode of the instant flow.

    % select a time point to present the instant flow.
    nt0 = 360; %360; %100; 500; 1000; ==>这些时刻的瞬时流场都试过了，但是统统不好看，索性就选第一个时刻（nt0=360）放在论文里面展示。
    %=====================================================================================================
    % Get flow field.
    load NearFlow_CF_sp8_NoBound_Tfine.mat;
    xg=linspace(xl_field(1),xl_field(2),Nx_field);
    yg=linspace(yl_field(1),yl_field(2),Ny_field); 
    zg=linspace(zl_field(1),zl_field(2),Nz_field); 
    [Xg,Yg,Zg]=meshgrid(xg,yg,zg);   
   
    if op==1 % raw flow.
        Ug0 = Ug(:,nt0);
        Vg0 = Vg(:,nt0);
        Wg0 = Wg(:,nt0);
        Ug0 = reshape(Ug0,size(Xg,1),size(Xg,2),size(Xg,3));
        Vg0 = reshape(Vg0,size(Xg,1),size(Xg,2),size(Xg,3));
        Wg0 = reshape(Wg0,size(Xg,1),size(Xg,2),size(Xg,3));
    elseif op==2 || op==3 % PCA-reconstructed flow.
        load NearFlow_CF_sp8_NoBound_Tfine_PCA_Accuracy99.99.mat;
        if op==2   % PCA reconstructed flow. 
            nt0_temp = 1+round((nt0-1)/2);
            Ug_PCA_rec = PCA_mode_Ug(:,1);
            Vg_PCA_rec = PCA_mode_Vg(:,1);
            Wg_PCA_rec = PCA_mode_Wg(:,1); 
            for i_nPCA=1:nmode  %15%
                Ug_PCA_rec = Ug_PCA_rec + Time_coef(nt0_temp,i_nPCA)*PCA_mode_Ug(:,i_nPCA+1);
                Vg_PCA_rec = Vg_PCA_rec + Time_coef(nt0_temp,i_nPCA)*PCA_mode_Vg(:,i_nPCA+1);
                Wg_PCA_rec = Wg_PCA_rec + Time_coef(nt0_temp,i_nPCA)*PCA_mode_Wg(:,i_nPCA+1);
            end 
        elseif op==3 % PCA mode. 
            i_mode = 4;  
            Ug_PCA_rec = PCA_mode_Ug(:,i_mode+1);  
            Vg_PCA_rec = PCA_mode_Vg(:,i_mode+1);  
            Wg_PCA_rec = PCA_mode_Wg(:,i_mode+1);                    
        end
        Ug0 = reshape(Ug_PCA_rec,size(Xg,1),size(Xg,2),size(Xg,3));
        Vg0 = reshape(Vg_PCA_rec,size(Xg,1),size(Xg,2),size(Xg,3)); 
        Wg0 = reshape(Wg_PCA_rec,size(Xg,1),size(Xg,2),size(Xg,3));                  
    end
               
    vel_raw=(Ug0.*Ug0+Vg0.*Vg0+Wg0.*Wg0).^0.5; 
    vel_nd=vel_raw.*2*pi;  %in units of 'flagellum length/ beat cycle',
    vel_dim = vel_nd.*scale_vel;  % in units of '\mu m/s'.
   
    
    
    %=====================================================================================================
    % Get sperm shape.
    load XNodes_sp8_NoBoundStokeslet_Tfine.mat;
    Q=Nhh+Ns;   
    nt_XNodes = size(Xs,2);
    Xtail_raw = Xs(Nhh+1:Q,:);  
    Ytail_raw = Xs(Q+Nhh+1:2*Q,:); 
    Ztail_raw = Xs(2*Q+Nhh+1:3*Q,:); 
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
    nt_flow = find(t_nd==t_nd_flow(end));
    Xtail_0 =  Xtail_sm(:,1:nt_flow);
    Ytail_0 =  Ytail_sm(:,1:nt_flow);
    Ztail_0 =  Ztail_sm(:,1:nt_flow);  
    Xtail_1=Xtail_0-repmat(Xtail_0(1,:),Ns,1);
    Ytail_1=Ytail_0-repmat(Ytail_0(1,:),Ns,1);
    Ztail_1=Ztail_0-repmat(Ztail_0(1,:),Ns,1);   
    Xtail_mean= mean(Xtail_1,2); %Ns*1  %The mean flagellar shape relative to the lab frame.
    Ytail_mean= mean(Ytail_1,2);
    Ztail_mean= mean(Ztail_1,2);
    dir_new = [1, 0, 0];
    dir_raw=[Xtail_mean(1)-Xtail_mean(end),Ytail_mean(1)-Ytail_mean(end),Ztail_mean(1)-Ztail_mean(end)];
    
    if op==1 || op==2 %instant sperm shape.
        nt0 = find(t_nd==t_nd_flow(nt0));
        % Aligned tail.
        Xtail_nt0 = Xtail_1(:,nt0); 
        Ytail_nt0 = Ytail_1(:,nt0); 
        Ztail_nt0 = Ztail_1(:,nt0); 
        [Xtail_align,Ytail_align,Ztail_align] = rxm_dir_align(dir_raw,dir_new,Xtail_nt0,Ytail_nt0,Ztail_nt0); %Xtail_align: Ns*1.
        % Aligned head.       
        head_tangent_nt0=head_tangent(nt0,:)';
        head_normal_nt0=head_normal(nt0,:)';
        head_binormal_nt0=cross(head_tangent_nt0,head_normal_nt0);
        B_nt0=[head_tangent_nt0,head_normal_nt0,head_binormal_nt0];
        [Bx_nt0_align,By_nt0_align,Bz_nt0_align] = rxm_dir_align(dir_raw,dir_new,B_nt0(1,:)',B_nt0(2,:)',B_nt0(3,:)');
        B_nt0_aligned = [Bx_nt0_align'; By_nt0_align'; Bz_nt0_align']; %3*3.  
        [xup1,xup2,xup3,xdown1,xdown2,xdown3,Ind_xup_neck] = Get_CF_Aligned_Head(B_nt0_aligned);  % xup1: M1*M2.   
    elseif op==3  %average sperm shape.
        % Aligned tail.
        [Xtail_align,Ytail_align,Ztail_align] = rxm_dir_align(dir_raw,dir_new,Xtail_mean,Ytail_mean,Ztail_mean); % Xtail_align: Ns*1. 
        % Aligned head.       
        B_mean=[1 0 0; 0 1 0; 0 0 1]; 
        [xup1,xup2,xup3,xdown1,xdown2,xdown3,Ind_xup_neck] = Get_CF_Aligned_Head(B_mean);  % xup1: M1*M2.   
    end   
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
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////
    view([170 15]); xlim([-1.5 0.5]); ylim([-1 1]); zlim([-1 1]); %view([190 10]); xlim([-1.5 0.5]); ylim([-0.8 0.8]); zlim([-0.8 0.8]); % for raw and PCA-recon. flow: nt=360.
    %
    %view([160 50]); %xlim([-1.5 0.5]); ylim([-1 1]); zlim([-0.4 0.6]); % for PCA mode 1.
    %view([160 20]); xlim([-1.5 0.5]); ylim([-1 1]); zlim([-1 1]); % for PCA mode 2.
    %view([160 20]); %xlim([-1.5 0.5]); ylim([-1 1]); zlim([-0.5 1]); % for PCA mode 3.
    %view([170 15]); %xlim([-1.5 0.5]); ylim([-1 1]); zlim([-0.8 0.8]); % for PCA mode 4.
    %view([180 20]); xlim([-1.5 0.5]); ylim([-1 1]); zlim([-1 1]); % for PCA mode 5.
    %view([156 14]); xlim([-1.5 0.5]); ylim([-1 1]); zlim([-1 1]); % for PCA mode 6.
    %xlim([-1.5 0.5]); ylim([-1 1]); zlim([-1 1]);  % for XYZ axes.
    
    
 if 1
    % Plot sperm.
    plot3(Xtail,Ytail,Ztail,'color',[36 100 171]./255,'linewidth',2);    %'color',[0.8500    0.3250    0.0980]
    CO_up(:,:,1) = 0.6*ones(size(xup1)); CO_up(:,:,2) = 0*ones(size(xup1)); CO_up(:,:,3) = 0*ones(size(xup1));
    CO_down(:,:,1) = 0*ones(size(xdown1)); CO_down(:,:,2) = 0.6*ones(size(xdown1)); CO_down(:,:,3) = 0*ones(size(xdown1));
    s1= surf(xup1,xup2,xup3,CO_up);
    s2= surf(xdown1,xdown2,xdown3,CO_down);
    s1.EdgeColor='none' ;
    s2.EdgeColor='none' ;
    
    % Plot flow streamlines.
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    [sx_3D,sy_3D,sz_3D] = meshgrid(xg(8:16),yg(7:12),zg(9:13)); % for raw and PCA-recon. flows.
    %[sx_3D,sy_3D,sz_3D] = meshgrid(xg(8:15),yg(1:21),zg(10:12)); % for PCA modes 1-2. 
    %[sx_3D,sy_3D,sz_3D] = meshgrid(xg(8:17),yg(6:17),zg(9:13)); % for PCA modes 3. 
    %[sx_3D,sy_3D,sz_3D] = meshgrid(xg(8:17),yg(7:14),zg(9:13)); % for PCA modes 4-6.  
    Generate_FlowStreamlines(Xg,Yg,Zg,Ug0,Vg0,Wg0,sx_3D,sy_3D,sz_3D, vel_dim);
    shading interp
    colorbar; 
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    caxis([0 20]);  % for raw and PCA-recon. flows. 
    %caxis([0 50]);  % for PCA mode 1. 
    %caxis([0 30]);  % for PCA modes 2-3. 
    %caxis([0 40]);  % for PCA modes 4-6. 
    
    camlight right
    lighting gouraud    
    return
 end    
 
    % Plot Cartesian axes.
    O = [-1 0 0]; % for PCA modes.
    Cx=[1; 0; 0]; 
    Cy=[0; 1; 0]; 
    Cz=[0; 0; 1]; 
    C=[Cx Cy Cz];
    Plot_Cartesian_Axes(fig,O,C)  
    camlight right
    lighting gouraud    
     
        
    return
end






%% Fig. 2: PCA cumulative variance.

if 0
    load FarFlow_CF_sp8_NoBound_Tfine_PCA_Accuracy99.99.mat;
    ind_90 = find(cumulative_variance>0.9,1,'first');
    %nmode = ind_90
    
    fs = 75;
    figure(2)
    clf; hold on; set(gcf,'color','w'); 
    axis square; box on;
    xlim([1 nmode]); ylim([0.2 1]);
    xticks(5:5:25);  yticks(0.2:0.2:1) %xticks([100:100:500]);
    %xticklabels({[],'20','40','60','80'}); 
    %text(7,0.23,num2str(ind_90),'FontName','times','FontSize',fs);
    %xticks(0:5:15); yticks(0.2:0.2:1)
    xlabel('PCA mode','FontName','times','FontSize',fs);  %,'interpreter','latex'
    ylabel('Cumulative variance','FontName','times','FontSize',fs);
    ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter = 'latex';
    
    
    plot(1:nmode, cumulative_variance(1:nmode), 'k','linewidth',5 )
    scatter(1:nmode, cumulative_variance(1:nmode), 300, 'k', 'filled'  )
    %plot([1 nmode],[0.9 0.9],'k:','linewidth',5)
    %plot([ind_90 ind_90],[0 1],'k:','linewidth',5)
    
    return
end





%% Fig. 3: time coefficients.

if 1

load NearFlow_CF_sp8_NoBound_Tfine.mat;
[~,np_bc] = min(abs(t_nd_flow-2*pi)); %np_bc~1000
% The time discretization in the flow_PCA data below is coarser than the
% flow data due to the computational cost limit, so we need to process
% their different 'np_bc'.
np_bc = round(np_bc/2);
np_rev = round((Freq_WF/Freq_rev)*np_bc);
n_rev = 8; 
n_bc_per_rev = round(np_rev/np_bc);  %n_bc_per_rev = round((Freq_WF/Freq_rev));


load NearFlow_CF_sp8_NoBound_Tfine_PCA_Accuracy99.99.mat;
nt = size(Time_coef,1); 
nmode = size(Time_coef,2);           
                
mode1 = 1;
mode2 = 2;
mode3 = 3;
mode4 = 4;

if 1 % mode1 versus mode2.
    fs = 85;
    figure(3)
    clf; set(gcf,'color','w');
    hold on; axis equal; box on;
    xlim([-0.37 0.65]);  ylim([-0.23 0.3]);  
    xticks(-0.2:0.2:0.6);  yticks(-0.2:0.2:0.2); 
    xlabel('1st PCA coefficient'); ylabel('2nd PCA coefficient') 
    set(gca,'TickLabelInterpreter','latex','FontName','times','FontSize',fs,'linewidth',1.5); 
    for i_rev = 1:n_rev % coefficient-plot envelope, over full revolution period.
        if i_rev==n_rev
            Time_coef_plot=Time_coef(1+(i_rev-1)*np_rev:i_rev*np_rev,:);
        else
            Time_coef_plot=Time_coef(1+(i_rev-1)*np_rev:2+i_rev*np_rev,:);
        end
        Time_coef_plot(end,2)=NaN;
        clrs=1:size(Time_coef_plot,1);
        p = patch(Time_coef_plot(:,mode1),Time_coef_plot(:,mode2),clrs,'EdgeColor','interp','linewidth',6);
        p.EdgeAlpha=0.3;
    end 
    for i_rev = 1:n_rev % over full revolution period.
        if i_rev==1
            Time_coef_plot=Time_coef(1+(i_rev-1)*np_rev:2+i_rev*np_rev,:);
            Time_coef_plot(end,2)=NaN;
            clrs=1:size(Time_coef_plot,1);
            p = patch(Time_coef_plot(:,mode1),Time_coef_plot(:,mode2),clrs,'EdgeColor','interp','linewidth',10);
            p.EdgeAlpha=1;
            for i_bc = 1:n_bc_per_rev+1  % over one revolution period.
                Time_coef_plot_temp=Time_coef(1+(i_bc-1)*np_bc:i_bc*np_bc,:);
                scatter(Time_coef_plot_temp(1,mode1),Time_coef_plot_temp(1,mode2),200*i_bc,[0.9 0 0],'filled')
            end
        end
        Time_coef_plot=Time_coef(1+(i_rev-1)*np_rev:i_rev*np_rev,:);        
        scatter(Time_coef_plot(np_rev,mode1),Time_coef_plot(np_rev,mode2),100*i_rev,[0.9 0 0],'filled','MarkerFaceAlpha',0.2)  
    end 
    return
    
else  % mode1 (model2) versus time.
    fs = 70;
    figure(3)
    clf; set(gcf,'color','w');
    hold on;  box on; daspect([1 0.0006 1])
    xlim([1 1*np_rev]);  ylim([-0.2 0.25]);   
    %xticks(1:np_rev:np_rev*n_rev+1); xticklabels({'0','1','2','3','4','5','6'}); yticks(-0.2:0.2:0.2); 
    %xlabel('t/T'); ylabel('Coefficients') 
    xtick_01 = np_bc*Freq_WF*0.1; xticks([1 xtick_01 2*xtick_01 3*xtick_01]);xticklabels({'0','0.1','0.2','0.3'});  %xticks(1:np_bc:n_bc_per_rev*np_bc+1); xticklabels({'0',[],[],num2str(3/Freq_WF)});  
    xlabel('Time (s)'); ylabel('Coefficients') 
    set(gca,'TickLabelInterpreter','latex','FontName','times','FontSize',1.2*fs,'linewidth',1.5);
    
    p1 = plot(1:10:nt,Time_coef(1:10:nt,mode1),'color',[36,100,171]/255,'linewidth',5); %'blue'
    p2 = plot(1:10:nt,Time_coef(1:10:nt,mode2),'color',[255 153 18]/255,'linewidth',5);    %'yellow'
    for i_bc = 1:n_bc_per_rev+1
        i_t = (i_bc-1)*np_bc+1;
        scatter(i_t,Time_coef(i_t,mode1),200*i_bc,[0.9 0 0],'filled') %'red'
        scatter(i_t,Time_coef(i_t,mode2),200*i_bc,[0.9 0 0],'filled') 
    end
    %scatter(1:np_bc:np_rev,Time_coef(1:np_bc:np_rev,mode1),600,[0.9 0 0],'filled')  %'red'
    %scatter(1:np_bc:np_rev,Time_coef(1:np_bc:np_rev,mode2),600,[0.9 0 0],'filled')
    lgd = legend([p1 p2],{'Mode 1','Mode 2'},'interpreter','latex'); 
    lgd.ItemTokenSize = [160 40];
    set(lgd,'FontSize',fs,'FontName','time','location','north','orientation','horizon');
    legend('boxoff')
    
end

return
end



%% Fig. 4: Colorbar.

figure(4)
clf; set(gcf,'color','w');
hold on; axis equal; axis off; box off;
    
colormap(parula); 
caxis([0,1]);

if 1
c = colorbar('Ticks',[0,1],'TickLabels',{'0','20'},'FontSize',48,'fontname','Times');  %,'FontSize',60  %,'TickLabels',{'0','1'}
c.Label.String = 'Flow velocity (\mum/s)';
c.TickLabelInterpreter='latex';
c.Position = [0.7 0.1 0.008 0.6];  % orientation: vertical
return
else
load NearFlow_CF_sp8_NoBound_Tfine.mat;
t_rev = 1/Freq_rev;
t_rev = roundn(t_rev,-2);
%c = colorbar('Ticks',[0,1],'TickLabels',{'0',num2str(t_rev)},'FontSize',60,'fontname','Times');  %,'TickLabels',{'0','1'}
c = colorbar('Ticks',[0,1],'TickLabels',{'0','0.30'},'FontSize',60,'fontname','Times');  %,'TickLabels',{'0','1'}
c.Label.String = 'Time (s)';
c.TickLabelInterpreter='latex';
c.Position = [0.85 0.2 0.01 0.7];  % orientation: vertical
end

c = colorbar('southoutside','Ticks',[0,1],'TickLabels',{'0','0.3'},'FontSize',80,'FontName','Times');  
c.Label.String = '(\mum/s)';
c.TickLabelInterpreter='latex'; 
c.Position = [0.2 0.5 0.45 0.015];  % orientation: horizontal

