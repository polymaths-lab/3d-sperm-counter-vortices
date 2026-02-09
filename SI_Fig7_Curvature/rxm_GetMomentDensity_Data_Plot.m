function rxm_GetMomentDensity_Data_Plot

% This function is to get the preliminary data for the power calculation
% (moment density) related to hydrodynamic force and internal stress (elasticity).
%
% Specifically, this function aims to get the following data along the flagellum:
% 1. bending and twisting of the flagellum along its local Frenet-Serret
% frames' vectors. --> for elastic moment density.
% 2. internal stress (integrated hydrodynamic force) along the flagellum
% local Frenet-Serret frames' vectors. --> for hydrodynamic moment density.
%
% Reference equations read from the paper: Rallabandi B, Wang Q, Potomkin M. Self-sustained three-dimensional beating of a model eukaryotic flagellum[J]. Soft Matter, 2022, 18(28): 5312-5322.


clc


%% Scaling parameters to dimensionalize the non-dimensional factors obtained from numerical reconstruction.


sp=8;
load 01_lab_frame_raw_traces_2017_2018_rxmSSppxyz.mat;
x_tail = X{sp};   
y_tail = Y{sp}; 
z_tail = Z{sp};    
x_tail = x_tail.*10^(-6);  %in units of m.
y_tail = y_tail.*10^(-6); 
z_tail = z_tail.*10^(-6); 

ns = size(x_tail,1);
nt = size(x_tail,2);
dt = 1/90;
t_exp = 0:dt:(nt-1)*dt; % in units of s.
a_temp = arclength(x_tail, y_tail, z_tail);  % size(a_temp)=ns*nt. in units of m.
s_mean = mean(a_temp,2); % size(s_mean)=ns*1. in units of m.
load FreeSperm_Frequency.mat;
w_fft = HF_Freq(sp);



visc = 0.001; %viscosity dynamics of the watery medium, in units of Pa*s.
L = s_mean(end); %characteristic length, equal to the measured human sperm's flagellar length, in units of m.
T = 1/w_fft; %characteristic beat period, equal to the measured human sperm's waveform beat cycle, in units of s.

scalep_u = L/T; %scaling parameter for velocity, in units of m/s.
scalep_f = visc*L/T; %scaling parameter for force density, in units of N/m.
scalep_m = visc*L*L/T; %scaling parameter for moment density, in units of N.
scalep_p = visc*L*L/T/T; %scaling parameter for power density, in units of W/m.
scalep_pp = 10^(-9); %scaling parameter to transfer the unit of power density from 'W/m' into 'fW/mum'.


clearvars -except scalep_u scalep_f  scalep_m  scalep_p scalep_pp L t_exp



%% Read smoothed sperm motiliy data.

load XNodes_sp8_NoBoundStokeslet.mat; %non-dimensional, both in time and space.
% xt123: lab-frame coordinates for sperm tail.
xt1 = xs(1+nhh:ns+nhh,:); %ns*nt. non-dimensional.
xt2 = xs(1+nhh+ns+nhh:ns+nhh+ns+nhh,:);
xt3 = xs(1+nhh+2*(ns+nhh):ns+nhh+2*(ns+nhh),:);
% extract 'ft123': force density along sperm tail, components along 3 lab-frame axes.
ft1 = f(1+nhh:ns+nhh,:);
ft2 = f(1+nhh+ns+nhh:ns+nhh+ns+nhh,:);
ft3 = f(1+nhh+2*(ns+nhh):ns+nhh+2*(ns+nhh),:);

nt = length(t_nd);
t_dim = t_nd.*t_exp(end)/t_nd(end);  %nt*1. in units of s.
dnt = round(length(t_nd)/ length(t_exp));
nBC = floor(t_nd(end)/(2*pi));
[~,nt_nBC] = min(abs(t_nd/(2*pi)-nBC)); % Locate the time point corresponding to the maximum beat cycle numbers.
s_mean = linspace(0,L,ns);



if 1 % Spatial smooth of the flagellar shape, to remove the noise arising from the raw experimental data being reconstructed.
    ppx = 0.999; % 'ppx=1' means interpolation without cubic spline.
    ppy = 0.999;
    ppz = 0.999;  
    for i_nt = 1:nt
        s_temp = arclength(xt1(:,i_nt), xt2(:,i_nt), xt3(:,i_nt));
        xt1(:,i_nt) = fnval(csaps(s_temp,xt1(:,i_nt),ppx),s_temp);
        xt2(:,i_nt) = fnval(csaps(s_temp,xt2(:,i_nt),ppy),s_temp);
        xt3(:,i_nt) = fnval(csaps(s_temp,xt3(:,i_nt),ppz),s_temp);
    end    
    xt1 = xt1-repmat(xt1(1,:),ns,1);
    xt2 = xt2-repmat(xt2(1,:),ns,1);
    xt3 = xt3-repmat(xt3(1,:),ns,1);    
end  
    
    

%% Get the flagellar bending moment density due to integrated hydrodynamic force.

ds = 1/(ns-1);  
ft1_int = nan(ns,nt);
ft2_int = nan(ns,nt);
ft3_int = nan(ns,nt);
ft1_int(ns-1,:) = (ft1(ns-1,:)+ft1(ns,:)).*ds/2;
ft2_int(ns-1,:) = (ft2(ns-1,:)+ft2(ns,:)).*ds/2;
ft3_int(ns-1,:) = (ft3(ns-1,:)+ft3(ns,:)).*ds/2;
for i_ns = 2:ns-1
    ft1_int(ns-i_ns,:) = ft1_int(ns-i_ns+1,:) + (ft1(ns-i_ns,:) + ft1(ns-i_ns+1,:)).*ds/2; %Non-NAN values range: (1:ns-1)*nt. non-dimensional.
    ft2_int(ns-i_ns,:) = ft2_int(ns-i_ns+1,:) + (ft2(ns-i_ns,:) + ft2(ns-i_ns+1,:)).*ds/2; 
    ft3_int(ns-i_ns,:) = ft3_int(ns-i_ns+1,:) + (ft3(ns-i_ns,:) + ft3(ns-i_ns+1,:)).*ds/2;       
end

ft_mag_int = (ft1_int.^2 + ft2_int.^2 + ft3_int.^2).^0.5; % magnitude of the total integrated hydrodynamic force (moment density).
ft_mag_int_dim = ft_mag_int.*scalep_m; %Non-NAN values range: (1:ns-1)*nt. in units of N.
 


%% Get the flagellar bending moment density due to internal elastic resistance.

kappa_mat = nan(size(xt1));
B0 = [0 0 1];
for nt0=1:nt
    X_temp_for_3D=[xt1(:,nt0), xt2(:,nt0), xt3(:,nt0)]; 
    [kappa_temp, ~, ~, ~, ~] = Get_Curvature_Torsion_NBT(X_temp_for_3D,B0);
            % [kappa, tau, tt, nn, bb] = Get_Curvature_Torsion_NBT_V2(x,B0)
            % kappa, tau: 1*ns array, but actually 'kappa' is 1*(ns-2) array and 'tau'
            % is 1*(ns-3) array because some values at the start and end points of the curve are 'NAN'.
            % kappa ... 1*(ns-2) vector, the measured curvature, with sign
            % tau ... 1*(ns-3) vector, the measured torsion, with sign
            % tt ... (ns-1)*3 array, measured tangent vectors, normalised
            % nn ... (ns-2)*3 array, measured normal vectors, normalised
            % bb ... (ns-2)*3 array, measured binormal vectors, normalised
     kappa_mat(:,nt0)= kappa_temp'; % Non-NAN values range: (2:ns-1)*nt. non-dimensional.     
end
kappa_mat_dim = kappa_mat./L; % Bending curvature, along the normal vector direction of the local frenet-serret frames: in units of '/m'.    
if 1 % Plot the curvature.
    fn='times';
    fs=62;
   
    figure(1) 
    clf; set(gcf,'color','w'); 
    hold on;  box on;
    set(gca,'position',[0.15 0.3 0.6 0.6]);  %daspect([0.8 60 1]);
    pcolor(t_dim(1:dnt:end),s_mean.*10^6,kappa_mat_dim(:,1:dnt:end)); 
    shading interp;
    kymo_temp_temp = kappa_mat_dim(:,1:dnt:end);      [min(kymo_temp_temp(:)),max(kymo_temp_temp(:))]
    xlim([0 t_dim(end)]); ylim([0 s_mean(end-1).*10^6]); 
    set(gca,'fontname',fn,'fontsize',fs);   %ax = gca;ax.FontSize = fs; ax.TickLabelInterpreter = 'latex'; 
    set(gca,'xtick',[0 1 2 3],'ytick',[0 20 40]); 
    hx=xlabel('Time (s)'); set(hx,'fontname',fn,'fontsize',fs+10);  %,'interpreter','latex'
    hy=ylabel('$s (\mu m)$'); set(hy,'interpreter','latex','fontname',fn,'fontsize',fs+10);  
    clrs0 = parula(2*64); clrs = clrs0(1:end,:);colormap(clrs);
    c=colorbar('AxisLocation','out'); c.Location='east'; c.FontSize=1*fs; c.FontName = 'times'; 
    c.Label.String = '\kappa (m^{-1})';  
    caxis([-3 3].*10^(5)); c.Ticks=[-3 3].*10^(5); c.TickLabels={'$-3\times10^{5}$','$3\times10^{5}$'}; 
    c.TickLabelInterpreter='latex';c.Label.FontSize=1*fs; 
    c.Position = [0.76 0.3 0.012 0.6]; 
    
    figure(11) 
    clf; set(gcf,'color','w'); 
    hold on;  box on;
    set(gca,'position',[0.15 0.3 0.6 0.6]);  %daspect([0.8 60 1]);
    pcolor(t_dim(1:dnt:end),s_mean.*10^6,abs(kappa_mat_dim(:,1:dnt:end))); 
    shading interp;
    xlim([0 t_dim(end)]); ylim([0 s_mean(end-1).*10^6]); 
    set(gca,'fontname',fn,'fontsize',fs);   %ax = gca;ax.FontSize = fs; ax.TickLabelInterpreter = 'latex'; 
    set(gca,'xtick',[0 1 2 3],'ytick',[0 20 40]); 
    hx=xlabel('Time (s)'); set(hx,'fontname',fn,'fontsize',fs+10);  %,'interpreter','latex'
    hy=ylabel('$s (\mu m)$'); set(hy,'interpreter','latex','fontname',fn,'fontsize',fs+10);  
    clrs0 = parula(2*64); clrs = clrs0(65:end,:);colormap(clrs);
    c=colorbar('AxisLocation','out'); c.Location='east'; c.FontSize=1*fs; c.FontName = 'times'; 
    c.Label.String = '|\kappa| (m^{-1})';  
    caxis([0 3].*10^(5)); c.Ticks=[0 3].*10^(5); c.TickLabels={'$0$','$3\times10^{5}$'}; 
    c.TickLabelInterpreter='latex';c.Label.FontSize=1*fs; 
    c.Position = [0.76 0.3 0.012 0.6]; 
    return
end



% Fit the bending stiffiness/ rigidity, 'E', for human sperm sperm.
% Reference paper:
% Gaffney E A, et al. Mammalian sperm motility: observation and theory, 2011. <-- 'E'-plot.
% Rallabandi B, et al. Self-sustained three-dimensional beating of a model eukaryotic flagellum, 2022. <-- anisotropic bending stiffiness values.
x_array = [0.25, 0.5, 1, 1.6, 2, 2.6, 3.5, 4.5, 5];
y_array = [7, 6.15, 5, 4, 3.65, 3, 2.45, 2, 1.85];
c = polyfit(x_array,y_array,5);
x = 0:0.01:5;
y = polyval(c,x);
E_temp = interp1(x,y,s_mean.*10^5,'spline'); 
E = E_temp.*10^(-21);  %1*ns. in units of Nm^2.
if 0% Check the polyfit curve.
    fn='times';
    fs=62;
    
    figure(2) 
    clf; set(gcf,'color','w','units','normalized') % make it white background
    hold on;axis equal;axis on;box on;axis square; %daspect([0.55 1 1]);  
    set(gca,'fontname','times','fontsize',fs); 
    %plot(x,y,'g','linewidth',6)
    %scatter(s_mean'.*10^5,E_temp,[],'r','filled')
    plot(s_mean'.*10^5,E_temp,'r','linewidth',6)
    xlim([0 5]); ylim([1 8])
    hx=xlabel('$s \ (10^{-5}m)$'); set(hx,'interpreter','latex','fontname',fn,'fontsize',fs+10);  
    hy=ylabel('$E \ (10^{-21}Nm^2)$'); set(hy,'interpreter','latex','fontname',fn,'fontsize',fs+10);  
    return
end


% Bending moment due to the passive elastic resistance, based on the curvature under the frenet-serret frames.
m_elas_bend = kappa_mat_dim(2:ns-1,:).*repmat(E(2:ns-1)',1,nt);  %(ns-2)*nt. in units of N*m.
f_elas_bend = Function_derivative_of_time(m_elas_bend',s_mean(2:ns-1),'time'); 
f_elas_bend = abs(f_elas_bend'); % moment density: (ns-2)*nt. in units of N.




%% Plot kymographs of the external and internal moment desities.

fn='times';
fs=68;

if 1
kymo_temp = ft_mag_int_dim;
figure(4) 
clf; set(gcf,'color','w'); 
hold on;  box on;
set(gca,'position',[0.15 0.3 0.5 0.5]); %set(gca,'position',[0.15 0.3 0.6 0.6]);  %daspect([0.8 60 1]);
pcolor(t_dim(1:dnt:end),s_mean(1:ns-1).*10^6,kymo_temp(1:ns-1,1:dnt:end)); 
shading interp;
kymo_temp_temp = kymo_temp(:,1:dnt:end);      [min(kymo_temp_temp(:)),max(kymo_temp_temp(:))]
xlim([0 t_dim(end)]); ylim([0 s_mean(end-1).*10^6]); xlim([0 1]);
set(gca,'fontname',fn,'fontsize',fs);   %ax = gca;ax.FontSize = fs; ax.TickLabelInterpreter = 'latex'; 
set(gca,'xtick',[0 1 2 3],'ytick',[0 20 40]); 
hx=xlabel('Time (s)'); set(hx,'fontname',fn,'fontsize',fs+10);  %,'interpreter','latex'
hy=ylabel('$s (\mu m)$'); set(hy,'interpreter','latex','fontname',fn,'fontsize',fs+10);  
clrs0 = parula(2*64); clrs = clrs0(64:end,:);colormap(clrs);
c=colorbar('AxisLocation','out'); c.Location='east'; c.FontSize=1*fs; c.FontName = 'times';
c.Label.String = '|m_{hyd}| (N)';  
caxis([0 2].*10^(-14)); c.Ticks=[0 2].*10^(-14); c.TickLabels={'$0$','$2\times10^{-14}$'}; 
c.TickLabelInterpreter='latex';c.Label.FontSize=1*fs; 
c.Position = [0.66 0.3 0.012 0.5];



kymo_temp = f_elas_bend;
figure(5)
clf; set(gcf,'color','w'); 
hold on;  box on;
set(gca,'position',[0.15 0.3 0.5 0.5]);   %set(gca,'position',[0.15 0.3 0.6 0.6]);  %daspect([0.8 60 1]);
pcolor(t_dim(1:dnt:end),s_mean(2:ns-1).*10^6,kymo_temp(:,1:dnt:end)); 
shading interp;
kymo_temp_temp = kymo_temp(:,1:dnt:end);      [min(kymo_temp_temp(:)),max(kymo_temp_temp(:))]
xlim([0 t_dim(end)]); ylim([0 s_mean(end-1).*10^6]);   xlim([0 1]);
set(gca,'fontname',fn,'fontsize',fs);   %ax = gca;ax.FontSize = fs; ax.TickLabelInterpreter = 'latex'; 
set(gca,'xtick',[0 1 2 3],'ytick',[0 20 40]); 
hx=xlabel('Time (s)'); set(hx,'fontname',fn,'fontsize',fs+10);  %,'interpreter','latex'
hy=ylabel('$s (\mu m)$'); set(hy,'interpreter','latex','fontname',fn,'fontsize',fs+10);  
clrs0 = parula(2*64); clrs = clrs0(65:end,:);colormap(clrs);
c=colorbar('AxisLocation','out'); c.Location='east'; c.FontSize=1*fs; c.FontName = 'times';
c.Label.String = '|m_{elas}| (N)';  
caxis([0 1].*10^(-10)); c.Ticks=[0 1].*10^(-10);c.TickLabels={'$0$','$1\times10^{-10}$'}; 
c.TickLabelInterpreter='latex';c.Label.FontSize=1*fs;
c.Position = [0.66 0.3 0.012 0.5];

return
end

%% Plot time average moment desities along the flagellum.

mean_f_temp = mean(ft_mag_int_dim(:,1:nt_nBC),2); max(mean_f_temp(:))
mean_f_temp = mean(f_elas_bend(:,1:nt_nBC),2); max(mean_f_temp(:))

if 0
figure(6) 
clf; set(gcf,'color','w'); 
hold on;  box on;
set(gca,'position',[0.3 0.3 0.5 0.31]);  %daspect([1.3 4*10^15 1]);
plot(s_mean.*10^6,mean(ft_mag_int_dim(:,1:nt_nBC),2),'linewidth',5)
xlim([0 s_mean(end-1).*10^6]); ylim([0 6].*10^(-15)); 
set(gca,'fontname',fn,'fontsize',fs);   %ax = gca;ax.FontSize = fs; ax.TickLabelInterpreter = 'latex'; 
set(gca,'ytick',[0 6].*10^(-15),'xtick',[0 20 40],'yticklabels',{'0','6\times 10^{-15}'}); 
hx=xlabel('$s (\mu m)$'); set(hx,'interpreter','latex','fontname',fn,'fontsize',fs+10);  
hy=ylabel('$\overline{|m_{hyd}|} (N)$'); set(hy,'interpreter','latex','fontname',fn,'fontsize',fs+10,'rotation',0);  



figure(7) 
clf; set(gcf,'color','w'); 
hold on;  box on;
set(gca,'position',[0.3 0.3 0.5 0.31]);   %set(gca,'position',[0.15 0.3 0.6 0.6]);  daspect([1.3 6*10^11 1]);
plot(s_mean(2:ns-1).*10^6,mean(f_elas_bend(:,1:nt_nBC),2),'linewidth',5)
xlim([0 s_mean(end-1).*10^6]);  ylim([0 4].*10^(-11)); 
set(gca,'fontname',fn,'fontsize',fs);   %ax = gca;ax.FontSize = fs; ax.TickLabelInterpreter = 'latex'; 
set(gca,'ytick',[0 4].*10^(-11),'xtick',[0 20 40],'yticklabels',{'0','4\times 10^{-11}'}); 
hx=xlabel('$s (\mu m)$'); set(hx,'interpreter','latex','fontname',fn,'fontsize',fs+10);  
hy=ylabel('$\overline{|m_{elas}|} (N)$'); set(hy,'interpreter','latex','fontname',fn,'fontsize',fs+10,'rotation',0);  


end