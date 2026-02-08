function rxm_AimFig_LabTrace_Plot


% This function plots the numerically reconstructed sperm trace (sp#8) in the lab frame. This is for the aim figure in our manuscript (MT_Fig1).
% The aim figure is for visualization purpose only, so the sperm body in this function is not to exact scale.

clc

%% Read sperm trace data.

% Experimental data.
sp=8;
load 01_lab_frame_raw_traces_2017_2018.mat; 
nt = size(Z{sp},2);
dt=1/90; % sampling frequency.
t_exp = 0:dt:(nt-1)*dt;
clearvars -except nt t_exp sp


% Numerical reconstruction data.
load XNodes_sp8_NoBoundStokeslet.mat;  %XNodes_sp8_BoundBlakelet_H0.2.mat; %
X = Xs(Nhh+1:Nhh+Ns,:); 
Y = Xs(2*Nhh+Ns+1:2*Nhh+2*Ns,:); 
Z = Xs(3*Nhh+2*Ns+1:3*Nhh+3*Ns,:);
t_dim = t_nd/(2*pi)/Freq_WF;
nt_fine=size(X,2); 
ns=size(X,1);
head_tangent=head_tangent';
head_normal=head_normal';
head_binormal=nan(size(head_tangent)); %3*nt_fine
for i_nt=1:nt_fine
    head_binormal(:,i_nt)=cross(head_tangent(:,i_nt),head_normal(:,i_nt));          
end


% Coarse-time interpolation, from numerical to experimental temporal
% points.
% Flagellum poition:
s_inp = linspace(0,1,ns);
[Tcoarse,Scoarse] = meshgrid(t_exp,s_inp); % query points for experimental coarse-time sperm tail coordinates.
[Tfine,Sfine] = meshgrid(t_dim,s_inp); % numerical fine-time sperm tail coordinates.
X = interp2(Tfine,Sfine,X,Tcoarse,Scoarse,'spline'); 
Y = interp2(Tfine,Sfine,Y,Tcoarse,Scoarse,'spline');
Z = interp2(Tfine,Sfine,Z,Tcoarse,Scoarse,'spline'); 
% Head orientation:
[Tcoarse,Scoarse] = meshgrid(t_exp,[1 2 3]) ;
[Tfine,Sfine] = meshgrid(t_dim,[1 2 3]) ; 
head_tangent = interp2(Tfine,Sfine,head_tangent,Tcoarse,Scoarse,'spline');             
head_normal  = interp2(Tfine,Sfine,head_normal,Tcoarse,Scoarse,'spline');             
head_binormal = interp2(Tfine,Sfine,head_binormal,Tcoarse,Scoarse,'spline');  


% Spatial smooth of the flagellar shape, to remove the noise arising from the raw experimental data being reconstructed.
ppx = 0.999; % 'ppx=1' means interpolation without cubic spline.
ppy = 0.999;
ppz = 0.999;  
for i_nt = 1:length(t_exp)
    s_temp = arclength(X(:,i_nt), Y(:,i_nt), Z(:,i_nt));
    X(:,i_nt) = fnval(csaps(s_temp,X(:,i_nt),ppx),s_temp);
    Y(:,i_nt) = fnval(csaps(s_temp,Y(:,i_nt),ppy),s_temp);
    Z(:,i_nt) = fnval(csaps(s_temp,Z(:,i_nt),ppz),s_temp);
end


clearvars -except nt ns X Y Z head_tangent head_normal head_binormal sp




%%  Generate sperm head (humanoid size).

    % creating mesh for sperm head provide by Hermes
    Ndiv = 60;
    radius_plot = 1.2/56;%normalized
    head_type = 1;    
    [Xhead,Yhead,Zhead] = sperm_head3D_shape(radius_plot, head_type, Ndiv);%size(Zhead)=[61 61]
    scal=1.4; %1.7;%scal=80;%scal=56;  %This 'scal' parameter can be adjusted according to visualization requirement.
    Xhead=Xhead*scal;Yhead=Yhead*scal;Zhead=Zhead*scal;
    M1=size(Xhead,1);M2=size(Xhead,2);M=M1*M2;
    [~,Ind_Neck] = max(Xhead(:)); Xhead_temp = Xhead(:); Xhead_neck = Xhead_temp(Ind_Neck);
    % head upside
    [zup_row,zup_col]=find(Zhead>=0); %length(zup_row)=1891
    Zhead_up=nan(size(Zhead));Xhead_up=nan(size(Xhead));Yhead_up=nan(size(Yhead));
    for i = 1:length(zup_row)
        Zhead_up(zup_row(i),zup_col(i))=Zhead(zup_row(i),zup_col(i));
        Xhead_up(zup_row(i),zup_col(i))=Xhead(zup_row(i),zup_col(i));
        Yhead_up(zup_row(i),zup_col(i))=Yhead(zup_row(i),zup_col(i));
    end
    Xhead_up = Xhead_up-Xhead_neck;
    xup=[reshape(Xhead_up,M,1);reshape(Yhead_up,M,1);reshape(Zhead_up,M,1)]; 
    % head downside
    [zdown_row,zdown_col]=find(Zhead<=0);%length(zdown_row)=1891
    Zhead_down=nan(size(Zhead));Xhead_down=nan(size(Xhead));Yhead_down=nan(size(Yhead));
    for i = 1:length(zdown_row)
        Zhead_down(zdown_row(i),zdown_col(i))=Zhead(zdown_row(i),zdown_col(i));
        Xhead_down(zdown_row(i),zdown_col(i))=Xhead(zdown_row(i),zdown_col(i));
        Yhead_down(zdown_row(i),zdown_col(i))=Yhead(zdown_row(i),zdown_col(i));
    end
    Xhead_down = Xhead_down-Xhead_neck;
    xdown=[reshape(Xhead_down,M,1);reshape(Yhead_down,M,1);reshape(Zhead_down,M,1)]; 



%% Plot sperm's lab-frame trace.

fig = 2;
figure(fig)
clf; set(gcf,'color','w','units','normalized');%,'outerposition',[0 0 1 1]) % make it white background
hold on; axis equal;
box off; axis off;
view(85,10); 
    
    
    % Read head spinning frequency.
    load FreeSperm_Frequency.mat;
    % (time in seconds)  
    dt = 1/90; % Each frame is taken every (1/90) seconds => 1.00 s corresponds to 90 video data frames (framerate was 90fps).
    tRange = 0:dt:(nt-1)*dt;
    ncycles = tRange(end)*HeadSpin_Freq(sp);%number of revolution cycles in the lab frame, during the recording time
    nspin = round(nt/ncycles); %nspin=13; nt=272; ncycle_spin=20;
    
    % Transparency parameter.
    FAlpha_mat_temp = linspace(-1,0,100);
    FAlpha_mat = 10.^FAlpha_mat_temp;
    ind_alpha = [30 50 85 90 95 100];
    i_alpha_temp = 0;
    
    for i_nt=[10*nspin+9, 13*nspin+9, 15*nspin, 16*nspin+3, 18*nspin+6, 19*nspin+7]  % for sp#8            
             % head surface rotation and translation
             B=[head_tangent(:,i_nt),head_normal(:,i_nt),head_binormal(:,i_nt)]; %orientation matrix (normalized) for sperm head
             neck=[X(2,i_nt),Y(2,i_nt),Z(2,i_nt)];  %neck=[X(1,i_nt),Y(1,i_nt),Z(1,i_nt)];               
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
             
             
             i_alpha_temp = i_alpha_temp+1;
             i_alpha = ind_alpha(i_alpha_temp);
             FAlpha = FAlpha_mat(i_alpha);
             % plot sperm head.
             s1= surf(head_up_x,head_up_y,head_up_z);
             s2= surf(head_down_x,head_down_y,head_down_z);
             cc = nan(size(s1.CData,1),size(s1.CData,2),3); cc(:,:,1)=177/255; cc(:,:,2)=24/255; cc(:,:,3)=45/255;s1.CData=cc; 
             cc = nan(size(s1.CData,1),size(s1.CData,2),3); cc(:,:,1)=218/255; cc(:,:,2)=207/255; cc(:,:,3)=168/255;s2.CData=cc; 
             s1.FaceAlpha = FAlpha;
             s2.FaceAlpha = FAlpha;
             % head projection.
             ss1= surf(head_up_x,head_up_y,head_up_z.*0-0.3 );
             ss2= surf(head_down_x,head_down_y,head_down_z.*0-0.3);
             cc = nan(size(s1.CData,1),size(s1.CData,2),3); cc(:,:,1)=0; cc(:,:,2)=0; cc(:,:,3)=0;ss1.CData=cc;
             cc = nan(size(s1.CData,1),size(s1.CData,2),3); cc(:,:,1)=0; cc(:,:,2)=0; cc(:,:,3)=0;ss2.CData=cc;             
             ss1.FaceAlpha = FAlpha;  
             ss2.FaceAlpha = FAlpha;
             %
             % plot sperm tail.
             PlotTailCylinder_AimFig(X(:,i_nt),Y(:,i_nt),Z(:,i_nt),FAlpha,fig)
             %plot3(X(:,i_nt),Y(:,i_nt),Z(:,i_nt),'color','r','linewidth',3) ;
    end
 
     
l1=light;
l1.Color = [1 1 1];
l1.Style = 'local';  %'infinite';  %
l1.Position = [0.8 -0.5 0.05];
l2=light;
l2.Color = [1 1 1];
l2.Style = 'local';  %'infinite';  %
l2.Position = [0 0 0.05];
l3=light;
l3.Color = [1 1 1];
l3.Style = 'local';  %'infinite';  %
l3.Position = [0.4 -0.2 0];
lighting gouraud   %phong  %flat  %
material metal  %shiny     

