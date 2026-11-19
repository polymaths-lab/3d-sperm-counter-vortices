function Trace_AbsDev_Plot

clc

% Read human sperm trace data.
load XNodes_sp8_NoBoundStokeslet.mat;
    human_trace_nb_temp = head_x0(1:10:end,:);  %nt*3
    human_trace_nb(:,1:2) = human_trace_nb_temp(:,1:2);
    human_trace_nb(:,3) = human_trace_nb_temp(:,3)+0.2;
load XNodes_sp8_BoundBlakelet_H0.2.mat;
    human_trace_wb = head_x0(1:10:end,:);  %nt*3
human_trace_AbsDev = sqrt(sum((human_trace_wb-human_trace_nb).^2,2));
t_bc_human = t_nd(1:10:end)./(2*pi);

% Read idealised sperm models.    
load Trace_VarySperm_100BC.mat;
i_head = 1;  % normal morphology.
ideal_trace_nb = [X_nb(:,i_head),  Y_nb(:,i_head),  Z_nb(:,i_head)];  % no-boundary trace 
ideal_trace_wb = [X_wb(:,i_head),  Y_wb(:,i_head),  Z_wb(:,i_head)];  % with-boundary trace 
ideal_trace_AbsDev = sqrt(sum((ideal_trace_nb-ideal_trace_wb).^2,2)); %nt*1. straight-line distance
t_bc_ideal = tRange/(2*pi);   
max(ideal_trace_AbsDev(:))
%return

figure(2)
    clf; set(gcf,'color','w');
    hold on; axis on; axis square; box on;
    p1 = plot(t_bc_human,human_trace_AbsDev,'g','linewidth',2);   
    p2 = plot(t_bc_ideal,ideal_trace_AbsDev,'b','linewidth',1); 
    ylim([-0.05 0.3]);
    xticks([0:20:100]); yticks([0:0.1:0.3]);
    fs=38;
    ax = gca;
    ax.FontSize = fs; 
    ax.TickLabelInterpreter = 'latex';
    ax.LineWidth = 1.5;  
    xlabel('Time','interpreter','latex','FontSize',fs+10);
    ylabel('$\Delta D$','interpreter','latex','FontSize',fs+10); 
    lgd = legend([p1 p2],{'Human','Ideal'}); 
    lgd.ItemTokenSize = [80 40];
    set(lgd,'FontSize',fs,'FontName','time','interpreter','latex','position',[.165 .42 .4 .1]);  %[.165 .39 .4 .1]);
    legend('boxoff')

%% Inset 1: traces of human sperm.        
    ax1 = axes('Position',[.25 .52 .4 .4]);
    hold on; axis equal; axis off; view(3); 
    set(gcf,'color','w');
    XL = [-0.1,1.1];YL = [-0.7,0.2];ZL = [-0.1,0.3];
    xlim([XL(1),XL(2)]);ylim([YL(1),YL(2)]);zlim([ZL(1),ZL(2)]);
    [xgrid,ygrid,zgrid] = meshgrid(linspace(XL(1),XL(2),10),...
                    linspace(YL(1),YL(2),10),...
                    linspace(ZL(1),ZL(2),10));                  
if 0% translucent plane color: blue.
   CO(:,:,1) = ones(10)*.0; % red 
   CO(:,:,2) = ones(10)*.4470; % green 
   CO(:,:,3) = ones(10)*0.7410; % blue 
elseif 0% translucent plane color: honeydew.
   CO(:,:,1) = ones(10)*240/255; % red 
   CO(:,:,2) = ones(10)*255/255; % green 
   CO(:,:,3) = ones(10)*240/255; % blue 
elseif 1
   CO(:,:,1) = ones(10)*0.8; % red 
   CO(:,:,2) = ones(10)*0.8; % green 
   CO(:,:,3) = ones(10)*0.8; % blue 
end       
surf(xgrid(:,:,1),ygrid(:,:,1),zgrid(:,:,1), CO,'FaceAlpha',0.05)% z-plane
surf(reshape(xgrid(end,:,:),[10,10]),reshape(ygrid(end,:,:),[10,10]),reshape(zgrid(end,:,:),[10,10]), CO,'FaceAlpha',0.1)%y plane
surf(reshape(xgrid(:,end,:),[10,10]),reshape(ygrid(:,end,:),[10,10]),reshape(zgrid(:,end,:),[10,10]), CO,'FaceAlpha',0.2)%x plane
shading interp

    p3 = plot3(human_trace_nb(:,1),human_trace_nb(:,2),human_trace_nb(:,3),'r','linewidth',2);
    p4 = plot3(human_trace_wb(:,1),human_trace_wb(:,2),human_trace_wb(:,3),'k','linewidth',2);
    
    plot3(human_trace_nb(:,1)*0+XL(2),human_trace_nb(:,2),human_trace_nb(:,3),'r','linewidth',0.5)
    plot3(human_trace_nb(:,1),human_trace_nb(:,2)*0+YL(2),human_trace_nb(:,3),'r','linewidth',0.5)
    plot3(human_trace_nb(:,1),human_trace_nb(:,2),human_trace_nb(:,3)*0+ZL(1),'r','linewidth',0.5)
    
    plot3(human_trace_wb(:,1)*0+XL(2),human_trace_wb(:,2),human_trace_wb(:,3),'k','linewidth',0.5)
    plot3(human_trace_wb(:,1),human_trace_wb(:,2)*0+YL(2),human_trace_wb(:,3),'k','linewidth',0.5)
    plot3(human_trace_wb(:,1),human_trace_wb(:,2),human_trace_wb(:,3)*0+ZL(1),'k','linewidth',0.5)
    
    lgd_ = legend([p3 p4],{'No boundary','Near boundary'}); 
    lgd_.ItemTokenSize = [60 40];  %lgd_.ItemTokenSize = [80 40];
    set(lgd_,'FontSize',fs-10,'FontName','time','interpreter','latex','position',[.46 .82 .4 .1]);
    legend('boxoff')

    
    
%% Inset 2: traces of idealised sperm models.        
    ax2 = axes('Position',[.38 .12 .4 .3]);
    hold on; axis equal; axis off; view(3); 
    set(gcf,'color','w');
    XL = [-2.1,0.2];YL = [-0.2,0.3];ZL = [-0.1,0.3];
    xlim([XL(1),XL(2)]);ylim([YL(1),YL(2)]);zlim([ZL(1),ZL(2)]);
    [xgrid,ygrid,zgrid] = meshgrid(linspace(XL(1),XL(2),10),...
                    linspace(YL(1),YL(2),10),...
                    linspace(ZL(1),ZL(2),10));                  
if 0% translucent plane color: blue.
   CO(:,:,1) = ones(10)*.0; % red 
   CO(:,:,2) = ones(10)*.4470; % green 
   CO(:,:,3) = ones(10)*0.7410; % blue 
elseif 0% translucent plane color: honeydew.
   CO(:,:,1) = ones(10)*240/255; % red 
   CO(:,:,2) = ones(10)*255/255; % green 
   CO(:,:,3) = ones(10)*240/255; % blue 
elseif 1% grey
   CO(:,:,1) = ones(10)*0.8; % red 
   CO(:,:,2) = ones(10)*0.8; % green 
   CO(:,:,3) = ones(10)*0.8; % blue 
end       
surf(xgrid(:,:,1),ygrid(:,:,1),zgrid(:,:,1), CO,'FaceAlpha',0.0)% z-plane
surf(reshape(xgrid(end,:,:),[10,10]),reshape(ygrid(end,:,:),[10,10]),reshape(zgrid(end,:,:),[10,10]), CO,'FaceAlpha',0.1)%y plane
surf(reshape(xgrid(:,end,:),[10,10]),reshape(ygrid(:,end,:),[10,10]),reshape(zgrid(:,end,:),[10,10]), CO,'FaceAlpha',0.2)%x plane
shading interp

    plot3(ideal_trace_nb(:,1),ideal_trace_nb(:,2),ideal_trace_nb(:,3),'r','linewidth',2);
    plot3(ideal_trace_wb(:,1),ideal_trace_wb(:,2),ideal_trace_wb(:,3),'k','linewidth',2);
    
    plot3(ideal_trace_nb(:,1)*0+XL(2),ideal_trace_nb(:,2),ideal_trace_nb(:,3),'r','linewidth',0.5)
    plot3(ideal_trace_nb(:,1),ideal_trace_nb(:,2)*0+YL(2),ideal_trace_nb(:,3),'r','linewidth',0.5)
    plot3(ideal_trace_nb(:,1),ideal_trace_nb(:,2),ideal_trace_nb(:,3)*0+ZL(1),'r','linewidth',0.5)
    
    plot3(ideal_trace_wb(:,1)*0+XL(2),ideal_trace_wb(:,2),ideal_trace_wb(:,3),'k','linewidth',0.5)
    plot3(ideal_trace_wb(:,1),ideal_trace_wb(:,2)*0+YL(2),ideal_trace_wb(:,3),'k','linewidth',0.5)
    plot3(ideal_trace_wb(:,1),ideal_trace_wb(:,2),ideal_trace_wb(:,3)*0+ZL(1),'k','linewidth',0.5)
        
    