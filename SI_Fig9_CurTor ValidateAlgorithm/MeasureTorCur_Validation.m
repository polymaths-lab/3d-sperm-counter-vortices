function MeasureTorCur_Validation

% This function validates my algorithm measuring the 3D curvature and torsion of a space curve (e.g. flagellar shape).

% Original flagellar shape and the corresponding reconstructed curve using
% the measured curvature and torsion are compared to validate the accuracy
% of the algorithm.

% My algorithm, meanwhile, is compared with Benjamin Friedrich's, which was
% used in the paper "Reconstruction of the three-dimensional beat pattern
% underlying swimming behaviors of sperm", and the code can be found at
% https://de.mathworks.com/matlabcentral/fileexchange/47885-frenet_robust-zip
% The comparison shows that mine has higher measure accuracy.

clc


%% Scaling parameter (arc length).

sp=8;
load 01_lab_frame_raw_traces_2017_2018_rxmSSppxyz.mat;
x_tail = X{sp};   
y_tail = Y{sp}; 
z_tail = Z{sp};    
x_tail = x_tail.*10^(-6);  %in units of m.
y_tail = y_tail.*10^(-6); 
z_tail = z_tail.*10^(-6); 
s = arclength(x_tail, y_tail, z_tail);  % ns*nt. in units of m.
s_mean = mean(s(end,:))*10^6;  %mean(s,2); % ns*1. in units of um.

clearvars -except s_mean


%% Load flagellar data--spatial denser and smoother.

load XNodes_sp8_NoBoundStokeslet.mat; % non-dimensional, both in time and space.
xt1_coarse = xs(1+nhh:ns+nhh,:); % ns*nt. lab-frame coordinates for sperm tail.
xt2_coarse = xs(1+nhh+ns+nhh:ns+nhh+ns+nhh,:);
xt3_coarse = xs(1+nhh+2*(ns+nhh):ns+nhh+2*(ns+nhh),:);

% Spatial smooth of the flagellar shape, to remove the noise arising from the raw experimental data.
    nt = length(t_nd);
    ppx = 0.999; % 'ppx=1' means interpolation without cubic spline.
    ppy = 0.999;
    ppz = 0.999;  
    ns_coarse = size(xt1_coarse,1);
    ns_fine = 100;
    for i_nt = 1:nt
        s_coarse = arclength(xt1_coarse(:,i_nt), xt2_coarse(:,i_nt), xt3_coarse(:,i_nt));         
        s_fine = interp1(linspace(0,1,ns_coarse),s_coarse,linspace(0,1,ns_fine));
        xt1_temp = interp1(linspace(0,1,ns_coarse),xt1_coarse(:,i_nt),linspace(0,1,ns_fine));
        xt2_temp = interp1(linspace(0,1,ns_coarse),xt2_coarse(:,i_nt),linspace(0,1,ns_fine));
        xt3_temp = interp1(linspace(0,1,ns_coarse),xt3_coarse(:,i_nt),linspace(0,1,ns_fine));
        
        xt1_fine(:,i_nt) = fnval(csaps(s_fine,xt1_temp',ppx),s_fine);
        xt2_fine(:,i_nt) = fnval(csaps(s_fine,xt2_temp',ppy),s_fine);
        xt3_fine(:,i_nt) = fnval(csaps(s_fine,xt3_temp',ppz),s_fine);
    end    
    xt1 = xt1_fine-repmat(xt1_fine(1,:),ns_fine,1);
    xt2 = xt2_fine-repmat(xt2_fine(1,:),ns_fine,1);
    xt3 = xt3_fine-repmat(xt3_fine(1,:),ns_fine,1);    
 
xt1 = s_mean*xt1(:,1:10:end); % ns*nt. dimensional lab-frame coordinates, in units of m.
xt2 = s_mean*xt2(:,1:10:end);
xt3 = s_mean*xt3(:,1:10:end);
nt = size(xt1,2);

%% Measure: curature, torsion, Frenet frames (T, N, B vectors).
   
nt0 = 10; %28;  % ==> for nt=310
%nt0 = 280; %99; ==> for nt=3100
Xcurve=[xt1(:,nt0), xt2(:,nt0), xt3(:,nt0)]; 


B0 = [0 0 1];
[kappa_rxm, tau_rxm, T_rxm, N_rxm, B_rxm] = Get_CurTor_TNB(Xcurve,B0); % kappa, tau: in units of '/m'.    


lwin = 4;  %lwin = 12; % moving span
weight=0; % regularization weight, default
%weight=0.1; % regularizes torsion (but introduces small artifacts)
[kappa_ben, tau_ben, T_ben, N_ben, B_ben]=frenet_robust(Xcurve',lwin,weight);
if 0 %nanmean(kappa_rxm'.*kappa_ben')<0
    kappa_ben=-kappa_ben;
end


%% Reconstruct the space curve using the measured curvature and torsion.


s_nt0 = arclength(Xcurve(:,1), Xcurve(:,2), Xcurve(:,3));  % ns*1. in units of m.

[~,~,~,~,Xcurve_rxm] = ReconCurve_by_CurTor(kappa_rxm(3:end-1), tau_rxm(3:end-1), T_rxm(3,:), N_rxm(3,:), Xcurve(3,:), s_nt0(3:end-1)); 
%[~,~,~,~,Xcurve_ben] = ReconCurve_by_CurTor(kappa_ben(8:end-6), tau_ben(8:end-6), T_ben(:,8)', N_ben(:,8)',Xcurve(8,:), s_nt0(8:end-6)); %lwin = 12;
[~,~,~,~,Xcurve_ben] = ReconCurve_by_CurTor(kappa_ben(4:end-1), tau_ben(4:end-1), T_ben(:,4)', N_ben(:,4)',Xcurve(4,:), s_nt0(4:end-1)); %lwin = 4;
    

%% Check measure accuracy: plot and deviation data.
if 0
    
fs=70;
fn='times';

figure(1) 
set(gcf,'color','w','units','normalized','outerposition',[0 0 1 1]) % make it white background
clf; hold on; axis equal;box on;
view([125 25]);  
ax = gca; ax.FontSize = fs; ax.FontName = fn; ax.TickLabelInterpreter = 'latex';

XL = [-55 5];YL = [-10,25];ZL = [-12 12];
xlim([XL(1),XL(2)]);ylim([YL(1),YL(2)]);zlim([ZL(1),ZL(2)]);
%xticks([-50,0]);yticks([-15,0,15]);zticks([-10,0,10])
xticks([]);yticks([]);zticks([])
   
ax = gca;ax.FontSize = 60; ax.TickLabelInterpreter = 'latex';
if 0
xlh=xlabel('$x (\mu m)$ ','interpreter','latex','FontSize',fs);
ylh=ylabel('$y (\mu m)$ ','interpreter','latex','FontSize',fs);
zlh=zlabel('$z (\mu m)$ ','interpreter','latex','FontSize',fs,'Rotation',0);
set(xlh,'position',[-25 YL(2) ZL(1)]);
set(ylh,'position',[XL(2)+10 10 ZL(1)-2]);
set(zlh,'position',[XL(2)+5 YL(1)-15 -5]);
end

%creating the side plane shading for 3D visualization======================
[xgrid,ygrid,zgrid] = meshgrid(linspace(XL(1),XL(2),10),...
                    linspace(YL(1),YL(2),10),...
                    linspace(ZL(1),ZL(2),10));                     
CO(:,:,1) = ones(10)*0;   %ones(10)*.0; % red
CO(:,:,2) = ones(10)*0;   %ones(10)*.4470; % green
CO(:,:,3) = ones(10)*0;   %ones(10)*0.7410; % blue
surf(xgrid(:,:,1),ygrid(:,:,1),zgrid(:,:,1), CO,'FaceAlpha',0.02)% z-plane
surf(reshape(xgrid(1,:,:),[10,10]),reshape(ygrid(1,:,:),[10,10]),reshape(zgrid(1,:,:),[10,10]), CO,'FaceAlpha',0.1)%y plane
surf(reshape(xgrid(:,1,:),[10,10]),reshape(ygrid(:,1,:),[10,10]),reshape(zgrid(:,1,:),[10,10]), CO,'FaceAlpha',0.2)%x plane
shading interp

% original flagellar shape.
p1 = plot3(Xcurve(:,1), Xcurve(:,2), Xcurve(:,3),'k','linewidth',5); %original curve
plot3(Xcurve(:,1)*0+XL(1),Xcurve(:,2),Xcurve(:,3),'k','linewidth',2); 
plot3(Xcurve(:,1),Xcurve(:,2)*0+YL(1),Xcurve(:,3),'k','linewidth',2); 
plot3(Xcurve(:,1),Xcurve(:,2),Xcurve(:,3)*0+ZL(1),'k','linewidth',2); 

if 1
% reconstructed flagellar shape (rxm).
p2 = plot3(Xcurve_rxm(:,1), Xcurve_rxm(:,2), Xcurve_rxm(:,3),'k:','linewidth',10); %reconstructed curve
plot3(Xcurve_rxm(:,1)*0+XL(1),Xcurve_rxm(:,2),Xcurve_rxm(:,3),'k:','linewidth',4);
plot3(Xcurve_rxm(:,1),Xcurve_rxm(:,2)*0+YL(1),Xcurve_rxm(:,3),'k:','linewidth',4);
plot3(Xcurve_rxm(:,1),Xcurve_rxm(:,2),Xcurve_rxm(:,3)*0+ZL(1),'k:','linewidth',4);

% local Frenet frames (rxm).
dns = 5;
quiver3(Xcurve(1:dns:end,1), Xcurve(1:dns:end,2), Xcurve(1:dns:end,3),T_rxm(1:dns:end,1),T_rxm(1:dns:end,2),T_rxm(1:dns:end,3),0.8,'g','linewidth',3,'MaxHeadSize',2)
quiver3(Xcurve(1:dns:end,1), Xcurve(1:dns:end,2), Xcurve(1:dns:end,3),N_rxm(1:dns:end,1),N_rxm(1:dns:end,2),N_rxm(1:dns:end,3),0.8,'b','linewidth',3,'MaxHeadSize',2)
quiver3(Xcurve(1:dns:end,1), Xcurve(1:dns:end,2), Xcurve(1:dns:end,3),B_rxm(1:dns:end,1),B_rxm(1:dns:end,2),B_rxm(1:dns:end,3),0.8,'r','linewidth',3,'MaxHeadSize',2)

else
    
% reconstructed flagellar shape (ben).
p2 = plot3(Xcurve_ben(:,1), Xcurve_ben(:,2), Xcurve_ben(:,3),'k:','linewidth',10); %reconstructed curve
plot3(Xcurve_ben(:,1)*0+XL(1),Xcurve_ben(:,2),Xcurve_ben(:,3),'k:','linewidth',4);
plot3(Xcurve_ben(:,1),Xcurve_ben(:,2)*0+YL(1),Xcurve_ben(:,3),'k:','linewidth',4);
plot3(Xcurve_ben(:,1),Xcurve_ben(:,2),Xcurve_ben(:,3)*0+ZL(1),'k:','linewidth',4);
 
% local Frenet frames (ben).
dns = 5;
quiver3(Xcurve(1:dns:end,1), Xcurve(1:dns:end,2), Xcurve(1:dns:end,3),T_ben(1,1:dns:end)',T_ben(2,1:dns:end)',T_ben(3,1:dns:end)',0.8,'g','linewidth',3,'MaxHeadSize',2)
quiver3(Xcurve(1:dns:end,1), Xcurve(1:dns:end,2), Xcurve(1:dns:end,3),N_ben(1,1:dns:end)',N_ben(2,1:dns:end)',N_ben(3,1:dns:end)',0.8,'b','linewidth',3,'MaxHeadSize',2)
quiver3(Xcurve(1:dns:end,1), Xcurve(1:dns:end,2), Xcurve(1:dns:end,3),B_ben(1,1:dns:end)',B_ben(2,1:dns:end)',B_ben(3,1:dns:end)',0.8,'r','linewidth',3,'MaxHeadSize',2)
  
end

if 0
legend(gca,[p1 p2],{'Original','Reconstructed'}...
    ,'FontName','times','FontSize',36,'ItemTokenSize',[120 40]);
legend('boxoff')
end



%return
end

%% Calculate measure accuracy/ deviation.
  
dev_rxm = nan(ns_fine,nt);
dev_ben = nan(ns_fine,nt);
i_dev = 1;
for i_nt = 1:1:nt
    Xcurve=[xt1(:,i_nt), xt2(:,i_nt), xt3(:,i_nt)]; 
    
    B0 = [0 0 1];
    [kappa_rxm, tau_rxm, T_rxm, N_rxm, B_rxm] = Get_CurTor_TNB(Xcurve,B0); % kappa, tau: in units of '/m'.    
    lwin = 4;  %lwin = 12; % moving span
    weight=0; % regularization weight, default
    %weight=0.1; % regularizes torsion (but introduces small artifacts)
    [kappa_ben, tau_ben, T_ben, N_ben, B_ben]=frenet_robust(Xcurve',lwin,weight);
    if 0 %nanmean(kappa_rxm'.*kappa_ben')<0
        kappa_ben=-kappa_ben;
    end
    
    s_nt0 = arclength(Xcurve(:,1), Xcurve(:,2), Xcurve(:,3));  % ns*1. in units of m.
    [~,~,~,~,Xcurve_rxm] = ReconCurve_by_CurTor(kappa_rxm(3:end-1), tau_rxm(3:end-1), T_rxm(3,:), N_rxm(3,:), Xcurve(3,:), s_nt0(3:end-1)); 
    %[~,~,~,~,Xcurve_ben] = ReconCurve_by_CurTor(kappa_ben(8:end-6), tau_ben(8:end-6), T_ben(:,8)', N_ben(:,8)',Xcurve(8,:), s_nt0(8:end-6)); %lwin = 12;
    [~,~,~,~,Xcurve_ben] = ReconCurve_by_CurTor(kappa_ben(4:end-1), tau_ben(4:end-1), T_ben(:,4)', N_ben(:,4)',Xcurve(4,:), s_nt0(4:end-1)); %lwin = 4;
    
    dis_rxm = Xcurve_rxm - Xcurve(3:end-1,:); % (ns-3)*3
    dev_rxm(3:end-1,i_dev) = sqrt(sum(dis_rxm.^2, 2)); % (ns-3)*1
    dis_ben = Xcurve_ben - Xcurve(4:end-1,:); % (ns-4)*3
    dev_ben(4:end-1,i_dev) = sqrt(sum(dis_ben.^2, 2)); % (ns-4)*1
    i_dev = i_dev+1;   
end

[max_rxm,ind] = max(mean(dev_rxm,2, 'omitnan')) % [2.2 99]
[max_ben,ind] = max(mean(dev_ben,2, 'omitnan')) % [3.7 99]

dev_rxm_mean = mean(dev_rxm(:), 'omitnan') %0.7
dev_ben_mean = mean(dev_ben(:), 'omitnan') %1.1

RMS_error_rxm = sqrt(mean(dev_rxm(~isnan(dev_rxm)).^2)) %1.45
RMS_error_ben = sqrt(mean(dev_ben(~isnan(dev_ben)).^2)) %2.20
% size(dev_rxm(~isnan(dev_rxm))) = [30070 1]

