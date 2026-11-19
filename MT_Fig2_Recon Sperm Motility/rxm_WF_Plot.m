function rxm_WF_Plot

% This function generates the observed flagellar waveform in the body frame of reference. 

clc


%% Load and process flagellum data.

sp=8;

load 03_Head_Fixed_Frame_Method_Rotations_Traces_raw_5_Microns_Difraction_head_spin_rxmSSppxyz.mat;
Xtail = X{sp}; 
Ytail = Y{sp}; 
Ztail = Z{sp}; 
nt=size(Xtail,2);
ns=size(Xtail,1);

for i_nt=1:nt
    arc_mat(:,i_nt) = arclength(Xtail(:,i_nt), Ytail(:,i_nt), Ztail(:,i_nt));
end
arc_mean = mean(arc_mat(end,:)) % measured as 46.35.

% Temporally interpolate the original waveform.
dt=1/90; %sampling frequency of the experimental imaging.
tRange=0:dt:(nt-1)*dt;
t0=linspace(0,tRange(end),nt);
t1=linspace(0,tRange(end),5*nt);
s=1:1:ns;
Xtail = interp2(t0,s',Xtail,t1,s'); 
Ytail = interp2(t0,s',Ytail,t1,s'); 
Ztail = interp2(t0,s',Ztail,t1,s'); 
nt=size(Xtail,2);
tRange=linspace(0,tRange(end),nt);

%% Generate sperm head: humanoid size.

    Ndiv = 60;
    radius_plot = 1.2/56; %normalized --> doesn't matter if head_type = 1.
    head_type = 1;    
    [Xhead,Yhead,Zhead] = sperm_head3D_shape(radius_plot, head_type, Ndiv);%size(Zhead)=[61 61]
    %Xrang=max(Xhead(:))-min(Xhead(:)); Yrang=max(Yhead(:))-min(Yhead(:)); Zrang=max(Zhead(:))-min(Zhead(:));
    %[Xrang/Zrang Yrang/Zrang Zrang/Zrang].*1.1 
    %return
    scal=1.2*arc_mean;  %In this plot, the head is for visualization, not used for reconstruction calculation, so the scale parameter is flexible.
    Xhead=Xhead*scal;
    Yhead=Yhead*scal;
    Zhead=Zhead*scal;
    Xneck=max(max(Xhead)); Xhead=Xhead-Xneck;
  
    
    M1=size(Xhead,1);M2=size(Xhead,2);M=M1*M2;
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
    
    
    
%% Plot the waveform figure.

figure(5) 
set(gcf,'color','w','units','normalized','outerposition',[0 0 1 1]) 
clf;hold on;
box on;axis equal;
view(-45,16); 

XL = [-7 60];YL = [-28 25];ZL = [-25 10];
xlim([XL(1),XL(2)]);ylim([YL(1),YL(2)]);zlim([ZL(1),ZL(2)]);
xticks([0,50]);yticks([-20,0,20]);zticks([-25,0])
   
fs=60;
ax = gca;ax.FontSize = fs; ax.TickLabelInterpreter = 'latex'; 
xlh=xlabel('$\xi_3 (\mu m)$ ','interpreter','latex','FontSize',fs);
ylh=ylabel('$\xi_1 (\mu m)$ ','interpreter','latex','FontSize',fs);
zlh=zlabel('$\xi_2 (\mu m)$ ','interpreter','latex','FontSize',fs,'Rotation',0);
set(xlh,'position',[10 YL(1)-2 ZL(1)]);
set(ylh,'position',[XL(1)-2 5 ZL(1)]-1);
set(zlh,'position',[XL(1)-2 YL(2)+7 -10]);
  

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
load FreeSperm_Frequency.mat;  % average beat frequency: 1/HF_Freq(sp)
ncycles = tRange(end)*HF_Freq(sp);             
nclrs = round(nt/ncycles)+1;  %nclrs = round(nt/ncycles); ==> '+1' is because only such can make the plotted period equal to the measured beat cycle determined by its frequency.
clrs = colormap(winter(nclrs)); 
nt_plot = 9*nclrs+1:10*nclrs;  %6*nclrs+1:7*nclrs
for i_nt = nt_plot
    if mod(i_nt,nclrs)==0
        plot3(Xtail(:,i_nt),Ytail(:,i_nt),Ztail(:,i_nt),'color',clrs(end,:),'linewidth',3) 
        plot3(Xtail(:,i_nt)*0+XL(2)-0.01,Ytail(:,i_nt),Ztail(:,i_nt),'color',clrs(end,:),'LineWidth',1); %tail: yz projection
    else
        plot3(Xtail(:,i_nt),Ytail(:,i_nt),Ztail(:,i_nt),'color',clrs(mod(i_nt,nclrs)+1,:),'linewidth',3) 
        plot3(Xtail(:,i_nt)*0+XL(2)-0.01,Ytail(:,i_nt),Ztail(:,i_nt),'color',clrs(mod(i_nt,nclrs)+1,:),'LineWidth',1); %tail: yz projection
    end
    plot3(Xtail(:,i_nt),Ytail(:,i_nt),Ztail(:,i_nt)*0+ZL(1)+0.01,'Color',[1 1 1]*0.5,'LineWidth',1); %tail: xy projection
    %plot3(Xtail(:,i_nt)*0+XL(2)-0.01,Ytail(:,i_nt),Ztail(:,i_nt),'Color',[1 1 1]*0.5,'LineWidth',1); %tail: yz projection
    plot3(Xtail(:,i_nt),Ytail(:,i_nt)*0+YL(2)-0.01,Ztail(:,i_nt),'Color',[1 1 1]*0.5,'LineWidth',1); %tail: xz projection
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
t_plot = tRange(nt_plot);
plot_period = t_plot(end)-t_plot(1);
plot_period = roundn(plot_period,-3);
colormap(winter); caxis([0,plot_period]);
c = colorbar('Ticks',[0,plot_period],'TickLabels',{'0',num2str(plot_period)},'FontSize',fs,'fontname','Times');  
c.Label.String = 'Time (s)'; c.TickLabelInterpreter='latex';
c.Position = [0.85 0.2 0.01 0.7]; 


