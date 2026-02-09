function rxm_IdealWF_Plot

% This function plots the flagellar waveform of the idealized sperm model. 

clc

%% Load and process the flagellum data.

[Xhead_up,Yhead_up,Zhead_up, Xhead_down,Yhead_down,Zhead_down] = Generate_Virtual_Head_UpDownSides;  %[M1 M2]=size(xup1);
xneck = max(Xhead_up(:));

load IdealWF_alpha0.0.mat;
Xtail = x(:,1:end-1);
Ytail = y(:,1:end-1);
Ztail = z(:,1:end-1);
nt = length(t)-1;
Xtail = Xtail+xneck; %size(x)=[ns nt]


    
%% Plot the waveform figure.

figure(1) 
set(gcf,'color','w','units','normalized','outerposition',[0 0 1 1]) 
clf;hold on;
box on;axis equal;
view(-65,15); 

XL = [-0.1,1.2];YL = [-0.3,0.3];ZL = [-0.2,0.2];
xlim([XL(1),XL(2)]);ylim([YL(1),YL(2)]);zlim([ZL(1),ZL(2)]);
xticks([0,1]);yticks([-0.3,0,0.3]);zticks([-0.2,0.2])
   
fs=90;%fs=70;
ax = gca;ax.FontSize = fs-20; ax.TickLabelInterpreter = 'latex'; 
xlh=xlabel('$\xi_3$ ','interpreter','latex','FontSize',fs);
ylh=ylabel('$\xi_1$ ','interpreter','latex','FontSize',fs);
zlh=zlabel('$\xi_2$ ','interpreter','latex','FontSize',fs,'Rotation',0);
%set(xlh,'position',[10 YL(1)-2 ZL(1)]);
%set(ylh,'position',[XL(1)-2 5 ZL(1)]-1);
%set(zlh,'position',[XL(1)-2 YL(2)+7 -10]);
  

[xgrid,ygrid,zgrid] = meshgrid(linspace(XL(1),XL(2),10),...
                    linspace(YL(1),YL(2),10),...
                    linspace(ZL(1),ZL(2),10));                     
CO(:,:,1) = ones(10)*1;   %ones(10)*.0; % red
CO(:,:,2) = ones(10)*1;   %ones(10)*.4470; % green
CO(:,:,3) = ones(10)*1;   %ones(10)*0.7410; % blue
surf(xgrid(:,:,1),ygrid(:,:,1),zgrid(:,:,1), CO,'FaceAlpha',0.05)% z-plane
surf(reshape(xgrid(end,:,:),[10,10]),reshape(ygrid(end,:,:),[10,10]),reshape(zgrid(end,:,:),[10,10]), CO,'FaceAlpha',0.15)%y plane
surf(reshape(xgrid(:,end,:),[10,10]),reshape(ygrid(:,end,:),[10,10]),reshape(zgrid(:,end,:),[10,10]), CO,'FaceAlpha',0.1)%x plane
shading interp
           

% Plot the head.
source_ligth = [90 80];
k = [0.9 1 1 5]; 
s1= surfl(Xhead_up,Yhead_up,Zhead_up,source_ligth,k );
s2= surfl(Xhead_down,Yhead_down,Zhead_down,source_ligth,k );
s1.EdgeColor='none' ;
s1.FaceColor=  [0.6 0 0 ];
s2.EdgeColor='none' ;
s2.FaceColor=  [0 0.6 0];
     


% Plot the tail.  
clrs = colormap(winter(nt)); 
for i_nt = 1:nt
    plot3(Xtail(:,i_nt),Ytail(:,i_nt),Ztail(:,i_nt),'color',clrs(i_nt,:),'linewidth',3) 
    plot3(Xtail(:,i_nt)*0+XL(2),Ytail(:,i_nt),Ztail(:,i_nt),'color',clrs(i_nt,:),'LineWidth',1); %tail: yz projection
    plot3(Xtail(:,i_nt),Ytail(:,i_nt),Ztail(:,i_nt)*0+ZL(1),'Color',[1 1 1]*0.5,'LineWidth',1); %tail: xy projection
    plot3(Xtail(:,i_nt),Ytail(:,i_nt)*0+YL(2),Ztail(:,i_nt),'Color',[1 1 1]*0.5,'LineWidth',1); %tail: xz projection
end


% Plot auxilary axes.
plot3(XL, [0 0], [0 0],'k','LineWidth',4); hold on;
plot3(XL, [0 0], [ZL(1) ZL(1)],'k','LineWidth',2);   %'color',[0.8500    0.3250    0.0980]
plot3([0 0], YL,[ZL(1) ZL(1)], 'k','LineWidth',2); 
plot3(XL, [YL(2) YL(2)], [0 0],'k','LineWidth',2); 
plot3([0 0], [YL(2) YL(2)], ZL,'k','LineWidth',2); 
plot3([XL(2) XL(2)], [0 0], ZL,'k','LineWidth',2); 
plot3([XL(2) XL(2)], YL, [0 0],'k','LineWidth',2); 


% Colorbar.
colormap(winter(500))
caxis([0 1])
c = colorbar('Ticks',[0,1],'TickLabels',{'0','1'},'FontSize',fs,'fontname','Times');  
c.Label.String = 'Time'; c.TickLabelInterpreter='latex';
c.Position = [0.85 0.2 0.02 0.7]; %c.Position = [0.85 0.2 0.01 0.7]; 


