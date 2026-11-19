function CurTor_Kymo_Plot
% plot measured curvature and torsion kymographs.

clc

load CurTor_sp8.mat;

   kappa_plot = kappa_rxm_mat;
   tau_plot = tau_rxm_mat;
   
   ns = size(tau_plot,1);
   nt = size(tau_plot,2);
   t_dim_mat = repmat(t_dim,ns,1);
   s_mean_mat = repmat(s_mean,1,nt);
   
   % select high-torsion-value samples.
   ind_tau_high = find(abs(tau_plot)>0.5);   % set the high-torsion-value limit as 0.1.
   ind_tau_high_selected = (t_dim_mat(ind_tau_high) >= 0.4) & (t_dim_mat(ind_tau_high) <= 0.6); % we only inspect a small range in 'time' to avoid messing the whole kymograph visualisation.
   ind_t_s = ind_tau_high(ind_tau_high_selected);
   
   figure(1)% tau 
   clf; set(gcf,'color','w','units','normalized','outerposition',[0 0 1 1]); 
   hold on; box on; %axis square;
   set(gca,'position',[0.15 0.3 0.5 0.5],'linewidth',1);  
   pcolor(t_dim,s_mean(3:end-1),tau_plot(3:end-1,:))
   scatter(t_dim_mat(ind_t_s),s_mean_mat(ind_t_s),50,'r','filled')
   shading interp;
    
   fn='times'; fs=68;
   xlim([0 1]); ylim([0 s_mean(end)]); %xlim([0 t_temp(end)]);
   set(gca,'fontname',fn,'fontsize',fs);   %ax = gca;ax.FontSize = fs; ax.TickLabelInterpreter = 'latex'; 
   set(gca,'xtick',[0 0.5 1],'ytick',[0 25 50]); 
   hx=xlabel('Time (s)'); set(hx,'interpreter','latex','fontname',fn,'fontsize',fs+10);  
   hy=ylabel('$s (\mu m)$'); set(hy,'interpreter','latex','fontname',fn,'fontsize',fs+10);  
   xticks([0 1]);yticks([0 20 40]);
    
   clrs0 = parula(2*64); clrs = clrs0(1:end,:);colormap(clrs);
   c=colorbar('AxisLocation','out'); c.Location='east'; c.FontSize=1*fs; c.FontName = 'times'; 
   %c.Label.String = '\tau (m^{-1})';  
   caxis([-0.3,0.3]); c.Ticks=[-0.3 0.3]; c.TickLabels={'-3','3'};  %c.TickLabels={'$-4\times10^{5}$','$4\times10^{5}$'}; 
   c.TickLabelInterpreter='latex';c.Label.FontSize=1*fs; 
   c.Position = [0.66 0.3 0.012 0.5];
    
   
 
   figure(2)% kappa 
   clf; set(gcf,'color','w','units','normalized','outerposition',[0 0 1 1]); 
   hold on; box on; %axis square;
   set(gca,'position',[0.15 0.3 0.5 0.5],'linewidth',1);  
   pcolor(t_dim,s_mean(2:end-1),kappa_plot(2:end-1,:))
   scatter(t_dim_mat(ind_t_s),s_mean_mat(ind_t_s),50,'r','filled')
   shading interp;
    
   fn='times'; fs=68;
   xlim([0 1]); ylim([0 s_mean(end)]); %xlim([0 t_temp(end)]);
   set(gca,'fontname',fn,'fontsize',fs);   %ax = gca;ax.FontSize = fs; ax.TickLabelInterpreter = 'latex'; 
   set(gca,'xtick',[0 0.5 1],'ytick',[0 25 50]); 
   hx=xlabel('Time (s)'); set(hx,'interpreter','latex','fontname',fn,'fontsize',fs+10);  
   hy=ylabel('$s (\mu m)$'); set(hy,'interpreter','latex','fontname',fn,'fontsize',fs+10);  
   xticks([0 1]);yticks([0 20 40]);
    
   clrs0 = parula(2*64); clrs = clrs0(1:end,:);colormap(clrs);
   c=colorbar('AxisLocation','out'); c.Location='east'; c.FontSize=1*fs; c.FontName = 'times'; 
   %c.Label.String = '\kappa (m^{-1})';  
   caxis([-0.3,0.3]); c.Ticks=[-0.3 0.3]; c.TickLabels={'-3','3'};  %c.TickLabels={'$-4\times10^{5}$','$4\times10^{5}$'}; 
   c.TickLabelInterpreter='latex';c.Label.FontSize=1*fs; 
   c.Position = [0.66 0.3 0.012 0.5];
      
   
    
   figure(3) % absolute kappa 
   clf; set(gcf,'color','w','units','normalized','outerposition',[0 0 1 1]); 
   hold on; box on; %axis square;
   set(gca,'position',[0.15 0.3 0.5 0.5],'linewidth',1);  
   pcolor(t_dim,s_mean(2:end-1),abs(kappa_plot(2:end-1,:)))
   scatter(t_dim_mat(ind_t_s),s_mean_mat(ind_t_s),50,'r','filled')
   shading interp;
    
   fn='times'; fs=68;
   xlim([0 1]); ylim([0 s_mean(end)]); %xlim([0 t_temp(end)]);
   set(gca,'fontname',fn,'fontsize',fs);   %ax = gca;ax.FontSize = fs; ax.TickLabelInterpreter = 'latex'; 
   set(gca,'xtick',[0 0.5 1],'ytick',[0 25 50]); 
   hx=xlabel('Time (s)'); set(hx,'interpreter','latex','fontname',fn,'fontsize',fs+10);  
   hy=ylabel('$s (\mu m)$'); set(hy,'interpreter','latex','fontname',fn,'fontsize',fs+10);  
   xticks([0 1]);yticks([0 20 40]);
    
   clrs0 = parula(2*64); clrs = clrs0(65:end,:);colormap(clrs);
   c=colorbar('AxisLocation','out'); c.Location='east'; c.FontSize=1*fs; c.FontName = 'times'; 
   %c.Label.String = '\kappa (m^{-1})';  
   caxis([0,0.3]); c.Ticks=[0 0.3]; c.TickLabels={'0','3'};  %c.TickLabels={'$-4\times10^{5}$','$4\times10^{5}$'}; 
   c.TickLabelInterpreter='latex';c.Label.FontSize=1*fs; 
   c.Position = [0.66 0.3 0.012 0.5];
      
   