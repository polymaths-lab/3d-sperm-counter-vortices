function Get_CurTor_Data


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
%s_mean = mean(s(end,:))*10^6;  %1*1. in units of um.
s_mean = mean(s,2)*10^6;  % ns*1. in units of um.

dt = 1/90;
nt = size(x_tail,2);
t_exp = 0:dt:(nt-1)*dt; % in units of s.

clearvars -except s_mean t_exp



%% Load flagellar data--spatial denser and smoother.

load XNodes_sp8_NoBoundStokeslet.mat; % non-dimensional, both in time and space.
xt1_coarse = xs(1+nhh:ns+nhh,:); % ns*nt. lab-frame coordinates for sperm tail.
xt2_coarse = xs(1+nhh+ns+nhh:ns+nhh+ns+nhh,:);
xt3_coarse = xs(1+nhh+2*(ns+nhh):ns+nhh+2*(ns+nhh),:);
t_dim = t_nd.*t_exp(end)/t_nd(end);  %nt*1. in units of s.

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
 
xt1 = s_mean(end)*xt1(:,1:10:end); % ns*nt. dimensional lab-frame coordinates, in units of m.
xt2 = s_mean(end)*xt2(:,1:10:end);
xt3 = s_mean(end)*xt3(:,1:10:end);
nt = size(xt1,2);
t_dim = t_dim(1:10:end);


%% Calculate curvature and torsion.
  

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
    
    kappa_rxm_mat(:,i_nt) = kappa_rxm(:);
    tau_rxm_mat(:,i_nt) = tau_rxm(:);
    kappa_ben_mat(:,i_nt) = kappa_ben(:);
    tau_ben_mat(:,i_nt) = tau_ben(:);  
end


X = [xt1; xt2; xt3];
save(['CurTor_sp8.mat'],'X','s_mean','t_dim','kappa_rxm_mat','tau_rxm_mat','kappa_ben_mat','tau_ben_mat');  
    