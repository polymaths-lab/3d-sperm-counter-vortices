function rxm_WF_video


clc

sp=8;


%% waveform interpolation.

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


X=Xtail; Y=Ytail; Z=Ztail;
X = X - repmat(X(1,:),ns,1);
Y = Y - repmat(Y(1,:),ns,1);
Z = Z - repmat(Z(1,:),ns,1);


%% generate sperm head: humanoid size

    Ndiv = 60;
    radius_plot = 1.2/56; %normalized --> doesn't matter if head_type = 1.
    head_type = 1;    
    [Xhead,Yhead,Zhead] = sperm_head3D_shape(radius_plot, head_type, Ndiv);%size(Zhead)=[61 61]
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
    xup=[reshape(Xhead_up,M,1);reshape(Yhead_up,M,1);reshape(Zhead_up,M,1)]; 
    % head downside
    [zdown_row,zdown_col]=find(Zhead<=0);%length(zdown_row)=1891
    Zhead_down=nan(size(Zhead));Xhead_down=nan(size(Xhead));Yhead_down=nan(size(Yhead));
    for i = 1:length(zdown_row)
        Zhead_down(zdown_row(i),zdown_col(i))=Zhead(zdown_row(i),zdown_col(i));
        Xhead_down(zdown_row(i),zdown_col(i))=Xhead(zdown_row(i),zdown_col(i));
        Yhead_down(zdown_row(i),zdown_col(i))=Yhead(zdown_row(i),zdown_col(i));
    end   
    xdown=[reshape(Xhead_down,M,1);reshape(Yhead_down,M,1);reshape(Zhead_down,M,1)]; 
    
    
%% waveform plot

% Font options for printing figures
fs=32;
fn='times';
fig = 4;

figure(fig) 
v=VideoWriter(['WF'],'MPEG-4');
v.FrameRate=30;
open(v);     

   
for i_nt=11:nt 
    
    set(gcf,'color','w','units','normalized','outerposition',[0 0 1 1]) % make it white background
    clf;hold on;
    box on;axis equal;
    view(-45,15);
    XL = [-10 60];YL = [-40 40];ZL = [-40 30]; % for sp8
    xlim([XL(1),XL(2)]);ylim([YL(1),YL(2)]);zlim([ZL(1),ZL(2)]);
    xticks([0:20:60]);yticks([-40:20:40]);zticks([-40:20:40])
    ax = gca;ax.FontSize = fs; ax.TickLabelInterpreter = 'latex';
    xlh=xlabel('$\xi_3 (\mu m)$','interpreter','latex','FontSize',1.2*fs);
    ylh=ylabel('$\xi_1 (\mu m)$','interpreter','latex','FontSize',1.2*fs);
    zlh=zlabel('$\xi_2 (\mu m)$','interpreter','latex','FontSize',1.2*fs,'Rotation',0);
    set(xlh,'position',[20 YL(1)-2 ZL(1)-6]);
    set(ylh,'position',[XL(1)-2 0 ZL(1)]-5);
    set(zlh,'position',[XL(1)-10 YL(2)+7 -10]);
    
    
    % plot sperm tail=========================================================
    if 0 %i_nt>11
        for k=1:i_nt-11
            PlotTailCylinder_Movie_LabTrajProj(X(:,k),Y(:,k),Z(:,k),0.1,fig, XL(2),YL(2),ZL(1))
        end  
    end
    for j = 1:nt %nclrs %
        pp = plot3(X(:,j),Y(:,j),Z(:,j)*0+ZL(1),'Color',[1 1 1]*0.85,'LineWidth',0.1); pp.Color(4)=0.3;
        pp = plot3(X(:,j)*0+XL(2),Y(:,j),Z(:,j),'Color',[1 1 1]*0.85,'LineWidth',0.1); pp.Color(4)=0.3;
        pp = plot3(X(:,j),Y(:,j)*0+YL(2),Z(:,j),'Color',[1 1 1]*0.85,'LineWidth',0.1); pp.Color(4)=0.3;
    end   
    % Transparency parameter.
    FAlpha_mat_temp = [-0.8:0.1:-0.4 -0.35:0.05:-0.15 0]; %FAlpha_mat_temp = linspace(-1,0,11);
    FAlpha_mat = 10.^FAlpha_mat_temp;               
    i_alpha=0;
    for k=i_nt-10:i_nt
       i_alpha = i_alpha+1;
       PlotTailCylinder_Movie_LabTraj(X(:,k),Y(:,k),Z(:,k),FAlpha_mat(i_alpha),fig, XL(2)-0.1,YL(2)-0.1,ZL(1)+0.1)
    end
    % plot tail mid-point trajectory and projection=========================================================
    if i_nt<21
        plot3(X(round(ns/2),1:i_nt),Y(round(ns/2),1:i_nt),Z(round(ns/2),1:i_nt),'color',[0.8500    0.3250    0.0980],'linewidth',1.5)   
        plot3(X(round(ns/2),1:i_nt)*0+XL(2)-0.1,Y(round(ns/2),1:i_nt),Z(round(ns/2),1:i_nt),'color',[0.8500    0.3250    0.0980],'linewidth',.1)
    else
        plot3(X(round(ns/2),i_nt-20:i_nt),Y(round(ns/2),i_nt-20:i_nt),Z(round(ns/2),i_nt-20:i_nt),'color',[0.8500    0.3250    0.0980],'linewidth',1.5)
        plot3(X(round(ns/2),i_nt-20:i_nt)*0+XL(2)-0.1,Y(round(ns/2),i_nt-20:i_nt),Z(round(ns/2),i_nt-20:i_nt),'color',[0.8500    0.3250    0.0980],'linewidth',.1)
    end
    scatter3(X(round(ns/2),i_nt),Y(round(ns/2),i_nt),Z(round(ns/2),i_nt),60,[0.8500    0.3250    0.0980],'filled')
    scatter3(X(round(ns/2),i_nt)*0+XL(2)-0.1,Y(round(ns/2),i_nt),Z(round(ns/2),i_nt),30,[0.8500    0.3250    0.0980],'filled')
      
    
    % Plot head.==================================================================
    % neck (head centre) position to locate/ translate the sperm head
    b1=[1;0;0];b2=[0;1;0];b3=[0;0;1];
    B=[b1,b2,b3];% B=[b1(:,i_nt),b2(:,i_nt),b3(:,i_nt)]; %orientation matrix (normalized) for sperm head
    neck=[X(1,i_nt),Y(1,i_nt),Z(1,i_nt)];          
    % head surface rotation and translation
    % upside head
    xupR=ApplyRotationMatrix(B,xup);
    xupT=TranslatePoints(xupR,neck);
    [head_up_x,head_up_y,head_up_z]=ExtractComponents(xupT);
    head_up_x=reshape(head_up_x,M1,M2);head_up_y=reshape(head_up_y,M1,M2);head_up_z=reshape(head_up_z,M1,M2);
    % downside head
    xdownR=ApplyRotationMatrix(B,xdown);
    xdownT=TranslatePoints(xdownR,neck);
    [head_down_x,head_down_y,head_down_z]=ExtractComponents(xdownT);
    head_down_x=reshape(head_down_x,M1,M2);head_down_y=reshape(head_down_y,M1,M2);head_down_z=reshape(head_down_z,M1,M2);
    % plot sperm head
    source_ligth = [90 80];
    k = [0.9 1 1 5]; % [ka kd ks shine]. By default, k is [.55 .6 .4 10].
    s1= surfl(head_up_x,head_up_y,head_up_z,source_ligth,k );
    s2= surfl(head_down_x,head_down_y,head_down_z,source_ligth,k );
    s1.EdgeColor='none' ;
    s1.FaceColor=  [177 24 45]/255;  %[0.6 0 0];  
    s2.EdgeColor='none' ;
    s2.FaceColor=  [218 207 168]/255;  %[0 0.6 0];
            
    
    % Plot auxilary axes.==================================================================
    plot3(XL, [0 0], [0 0],'k','LineWidth',2); hold on;
    plot3(XL, [0 0], [ZL(1) ZL(1)],'k','LineWidth',1);   
    plot3([0 0], YL,[ZL(1) ZL(1)], 'k','LineWidth',1); 
    plot3(XL, [YL(2) YL(2)], [0 0],'k','LineWidth',1); 
    plot3([0 0], [YL(2) YL(2)], ZL,'k','LineWidth',1); 
    plot3([XL(2) XL(2)], [0 0], ZL,'k','LineWidth',1); 
    plot3([XL(2) XL(2)], YL, [0 0],'k','LineWidth',1); 
   

    % Title: time.==================================================================      
    inst = dt/5*(i_nt-1);
    format short
    inst = roundn(inst,-2);
    title(['Time: ',num2str(inst),'s'],'Fontname','Times New Roman','fontsize',fs,'fontweight','normal');

    
    l1=light;
    l1.Color = [1 1 1];
    l1.Style = 'local';  %'infinite';  %
    l1.Position = [20 0 20];
    lighting gouraud   %phong  %flat  %
    material metal  %shiny     

        
    drawnow;
    frame(i_nt)=getframe(gcf);
    writeVideo(v,frame(i_nt));

end

end
