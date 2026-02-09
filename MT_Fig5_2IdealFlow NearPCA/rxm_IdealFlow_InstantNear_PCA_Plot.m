function rxm_IdealFlow_InstantNear_PCA_Plot

% This function plots the instant flow field, as well as its PCA reconstruction, around the reconstructed human sperm, relative to the comoving frame of reference.
% Figures:
% 1. 3D flow (raw and PCA-recon., and PCA modes).
% 2. PCA cumulative variance.
% 3. PCA time coefficients.
% 4. Colorbar.

clc


%% Fig. 1: 3D flow.

if 1
    op = 3;
    % op=1: plot the raw instant flow.
    % op=2: plot the PCA-reconstructed instant flow.
    % op=3: plot the PCA mode of the instant flow.


    % select a time point to present the instant flow.
    nt0 = 350; % for the dataset of 'NearFlow_CF_alpha0.6_Tfine.mat'
    %nt0 = 70; % for the dataset of 'NearFlow_CF_alpha0.6.mat'
    %=====================================================================================================
    % Get flow field.
    load NearFlow_CF_alpha0.6_Tfine.mat;
    %load NearFlow_CF_alpha0.6.mat;
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
        load NearFlow_CF_alpha0.6_Tfine_PCA_Accuracy99.99.mat;
        if op==2   % PCA reconstructed flow.  
            Ug_PCA_rec = PCA_mode_Ug(:,1);
            Vg_PCA_rec = PCA_mode_Vg(:,1);
            Wg_PCA_rec = PCA_mode_Wg(:,1); 
            for i_nPCA=1:7%nmode 
                Ug_PCA_rec = Ug_PCA_rec + Time_coef(nt0,i_nPCA)*PCA_mode_Ug(:,i_nPCA+1);
                Vg_PCA_rec = Vg_PCA_rec + Time_coef(nt0,i_nPCA)*PCA_mode_Vg(:,i_nPCA+1);
                Wg_PCA_rec = Wg_PCA_rec + Time_coef(nt0,i_nPCA)*PCA_mode_Wg(:,i_nPCA+1);
            end 
        elseif op==3 % PCA mode. 
            i_mode = 1; 
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
    
    
    %=====================================================================================================
    % Get sperm shape.
    load XNodes_Lab_alpha0.6_Tfine.mat;
    %load XNodes_Lab_alpha0.6.mat;
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

    if op==1 || op==2 % instant sperm shape.
        nt0_X = find(tRange==tRange_flow(nt0));
        % Aligned tail.
        Xtail_nt0 = Xtail_1(:,nt0_X); 
        Ytail_nt0 = Ytail_1(:,nt0_X); 
        Ztail_nt0 = Ztail_1(:,nt0_X); 
        [Xtail_align,Ytail_align,Ztail_align] = rxm_dir_align(dir_raw,dir_new,Xtail_nt0,Ytail_nt0,Ztail_nt0); %Xtail_align: Ns*1.
        % Aligned head.       
        head_tangent_nt0=head_tangent(nt0_X,:)';
        head_normal_nt0=head_normal(nt0_X,:)';
        head_binormal_nt0=cross(head_tangent_nt0,head_normal_nt0);
        B_nt0=[head_tangent_nt0,head_normal_nt0,head_binormal_nt0]; 
        [Bx_nt0_align,By_nt0_align,Bz_nt0_align] = rxm_dir_align(dir_raw,dir_new,B_nt0(1,:)',B_nt0(2,:)',B_nt0(3,:)');
        B_nt0_aligned = [Bx_nt0_align'; By_nt0_align'; Bz_nt0_align']; %3*3.  
        [xup1,xup2,xup3,xdown1,xdown2,xdown3,Ind_xup_neck] = Get_CF_Aligned_Head(B_nt0_aligned);  % xup1: M1*M2.          
    elseif op==3 % average sperm shape.
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
    fig = 4;
    
    figure(fig)
    clf; set(gcf,'color','w');
    hold on; axis equal; %axis off;     
    %/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    view(3); xlim([xl_field(1),xl_field(2)]); ylim([yl_field(1),yl_field(2)]); zlim([zl_field(1),zl_field(2)]); % for PCA modes 1-2.
    %view([-75 20]); xlim([xl_field(1),xl_field(2)]); ylim([yl_field(1),yl_field(2)]); zlim([zl_field(1),zl_field(2)]); % for PCA modes 3-4.
    %view([-20 25]); xlim([xl_field(1),xl_field(2)]); ylim([yl_field(1),yl_field(2)]); zlim([zl_field(1),zl_field(2)]); % for PCA modes 5-6.
    %view([-20 25]); xlim([xl_field(1),xl_field(2)]); ylim([yl_field(1),yl_field(2)]); zlim([zl_field(1),zl_field(2)]); % for raw flow and its recon: nt0=350.
     

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
    % for raw and PCA-recon. flows.
    [sx_3D,sy_3D,sz_3D] = meshgrid(xg(1:1:end),yg(7:2:15),zg(9:2:13)); % for PCA modes 1-2.
    %[sx_3D,sy_3D,sz_3D] = meshgrid(xg(1:1:end),yg(3:2:19),zg(2:9:20)); % for PCA modes 3-4.
    %[sx_3D,sy_3D,sz_3D] = meshgrid(xg(1:1:end),yg(9:13),zg(9:13)); % for PCA modes 5-6.
    %[sx_3D,sy_3D,sz_3D] = meshgrid(xg(1:1:end),yg(6:16),zg(10:12)); % for raw flow and its recon.: nt0=350.
    Generate_FlowStreamlines(Xg,Yg,Zg,Ug0,Vg0,Wg0,sx_3D,sy_3D,sz_3D, vel_nd);
    shading interp
    colorbar; 
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    caxis([0 0.1]);  % for PCA modes.
    %caxis([0 0.01]);  % for the raw and recon. flows: nt0=350.
        
    camlight right
    lighting gouraud    
    return
 end    
 
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




%% Fig. 2: PCA cumulative variance.

if 1
    load NearFlow_CF_alpha0.6_Tfine_PCA_Accuracy99.99.mat;
    %ind_90 = find(cumulative_variance>0.9,1,'first');
    nmode = ind_90;
    
    fs = 74;
    figure(1)
    clf; hold on; set(gcf,'color','w'); 
    axis square; box on;
    xlim([1 nmode]); ylim([0.3 1]);
    %xticks([ind_90 20:20:80]); yticks([0.4:0.2:0.8 0.9 1])
    %xticklabels({[],'20','40','60','80'}); 
    %text(7,0.23,num2str(ind_90),'FontName','times','FontSize',fs);
    xticks(1:7); yticks(0.4:0.2:1)
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

if 0

load NearFlow_CF_alpha0.6_Tfine_PCA_Accuracy99.99.mat;
np_rev = 600;
np_bc = 100;
n_rev = 16;
n_bc_per_rev = 6;
%load NearFlow_CF_alpha0.6_PCA_Accuracy99.99.mat;
%np_rev = 120;
%np_bc = 20;
%n_rev = 1;
%n_bc_per_rev = 6;
nt = size(Time_coef,1); 
nmode = size(Time_coef,2); 

mode1 = 1;
mode2 = 2;  

if 0 % mode1 versus mode2.
    fs = 68;
    figure(2)
    clf; set(gcf,'color','w');
    hold on; axis equal; box on;
    plot(Time_coef(:,mode1),Time_coef(:,mode2))
    return
    xlim([-0.11 0.11]);  ylim([-0.11 0.11]);  
    xticks(-0.1:0.1:0.1);  yticks(-0.1:0.1:0.1); 
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
        %p.EdgeAlpha=0.05;
        p.EdgeAlpha=1;
    end 
    if 0
    for i_rev = 1:n_rev % over full revolution period.
        Time_coef_plot=Time_coef(1+(i_rev-1)*np_rev:i_rev*np_rev,:);
        if i_rev==1
            Time_coef_plot=Time_coef(1+(i_rev-1)*np_rev:2+i_rev*np_rev,:);
            Time_coef_plot(end,2)=NaN;
            clrs=1:size(Time_coef_plot,1);
            p = patch(Time_coef_plot(:,mode1),Time_coef_plot(:,mode2),clrs,'EdgeColor','interp','linewidth',5);
            p.EdgeAlpha=1;
            for i_bc = 1:n_bc_per_rev  % over one revolution period.
                Time_coef_plot_temp=Time_coef(1+(i_bc-1)*np_bc:i_bc*np_bc,:);
                scatter(Time_coef_plot_temp(40,mode1),Time_coef_plot_temp(40,mode2),200*i_bc,[0.9 0 0],'filled')
            end
        end
        scatter(Time_coef_plot(572,mode1),Time_coef_plot(572,mode2),100*i_rev,[0.9 0 0],'filled','MarkerFaceAlpha',0.2)  
    end
    end
    return
    
else  % mode1 (model2) versus time.
    fs = 78;
    figure(3)
    clf; set(gcf,'color','w');
    hold on;  box on; 
    plot(1:1*np_rev,Time_coef(1:1*np_rev,mode1))
    plot(1:1*np_rev,Time_coef(1:1*np_rev,mode2))
    return
    
    daspect([1 0.0007 1])
    xlim([1 1*np_rev+1]);  ylim([-0.11 0.11]);  
    xticks(1:100:1*np_rev+1); xticklabels({'0','1','2','3','4','5','6'}); yticks(-0.1:0.1:0.1); 
    xlabel('t/T'); ylabel('Coefficients') 
    set(gca,'TickLabelInterpreter','latex','FontName','times','FontSize',fs,'linewidth',1.5);
    p1 = plot(1:1*np_rev,Time_coef(1:1*np_rev,mode1),'color',[36/255,100/255,171/255],'linewidth',5); %'blue'
    p2 = plot(1:1*np_rev,Time_coef(1:1*np_rev,mode2),'color',[255 153 18]/255,'linewidth',5);    %'yellow'
    for i_bc = 1:n_bc_per_rev%+1
        i_t = (i_bc-1)*np_bc+1;
        scatter(i_t,Time_coef(i_t,mode1),200*i_bc,[0.9 0 0],'filled') %'red'
        scatter(i_t,Time_coef(i_t,mode2),200*i_bc,[0.9 0 0],'filled') 
    end    
    %scatter(40:np_bc:np_rev,Time_coef(40:np_bc:np_rev,mode1),600,[0.9 0 0],'filled')
    %scatter(40:np_bc:np_rev,Time_coef(40:np_bc:np_rev,mode2),600,[0.9 0 0],'filled')
    lgd = legend([p1 p2],{'Mode 1','Mode 2'},'interpreter','latex'); 
    lgd.ItemTokenSize = [160 40];
    set(lgd,'FontSize',fs,'FontName','time','location','north','orientation','horizon');
    legend('boxoff')
end

return
end


%% Fig. 4: Colorbar.

figure(2)
clf; set(gcf,'color','w');
hold on; axis equal; axis off; box off;
    
colormap(parula); 
caxis([0,1]);

%c = colorbar('Ticks',[0,1],'TickLabels',{'0','0.01'},'FontSize',50,'fontname','Times');  %,'TickLabels',{'0','1'}
%c.Label.String = 'Flow velocity';
c = colorbar('Ticks',[0,1],'TickLabels',{'0','6'},'FontSize',50,'fontname','Times');  %,'TickLabels',{'0','1'}
c.Label.String = 'Time';
c.TickLabelInterpreter='latex';
c.Position = [0.85 0.2 0.008 0.55];  % orientation: vertical

return
c = colorbar('southoutside','Ticks',[0,1],'TickLabels',{'0','0.3'},'FontSize',80,'FontName','Times');  
c.Label.String = 'Flow velocity';
c.TickLabelInterpreter='latex'; 
c.Position = [0.2 0.3 0.45 0.015];  % orientation: horizontal


c = colorbar('southoutside','Ticks',[0,1],'TickLabels',{'0','0.3'},'FontSize',80,'FontName','Times');  
c.Label.String = '(\mum/s)';
c.TickLabelInterpreter='latex'; 
c.Position = [0.2 0.5 0.45 0.015];  % orientation: horizontal

