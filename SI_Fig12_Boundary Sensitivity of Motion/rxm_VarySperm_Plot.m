function rxm_VarySperm_Plot


clc

if 0
load Sperm_Head1_Alpha0.mat;
fig=1;

figure(fig)
clf; set(gcf,'color','w','units','normalized');
hold on; axis equal; box off; axis off;
view(3); 

% plot sperm head.
s1= surf(xup1,xup2,xup3);
s2= surf(xdown1,xdown2,xdown3);
cc = nan(size(s1.CData,1),size(s1.CData,2),3); cc(:,:,1)=177/255; cc(:,:,2)=24/255; cc(:,:,3)=45/255; s1.CData=cc; 
cc = nan(size(s2.CData,1),size(s2.CData,2),3); cc(:,:,1)=218/255; cc(:,:,2)=207/255; cc(:,:,3)=168/255; s2.CData=cc;           

% plot sperm tail.
% Transparency parameter.
FAlpha_mat_temp = linspace(-1,0,21);
FAlpha_mat = 100.^FAlpha_mat_temp;
i_alpha_temp = 0;
for i_t = 1:size(x,2)    
    i_alpha_temp = i_alpha_temp+1;
    FAlpha = FAlpha_mat(i_alpha_temp);
    p = plot3(x(:,i_t),y(:,i_t),z(:,i_t),'color',[36 100 171]./255,'linewidth',4);    %'color',[0.8500    0.3250    0.0980]
    p.Color(4) = FAlpha;
    %PlotTailCylinder(x(:,i_t),y(:,i_t),z(:,i_t),FAlpha,fig)  
end           
shading interp
    
l1=light;
l1.Color = [1 1 1];
l1.Style = 'local';  %'infinite';  %
l1.Position = [-0.1 -0.3 -0.1];
l2=light;
l2.Color = [1 1 1];
l2.Style = 'local';  %'infinite';  %
l2.Position = [0.4 -0.3 0.1];  
l3=light;
l3.Color = [1 1 1];
l3.Style = 'local';  %'infinite';  %
l3.Position = [0.6 0 0.3];
l4.Color = [1 1 1];
l4.Style = 'local';  %'infinite';  %
l4.Position = [0.8 0 0.3];
lighting gouraud   %phong  %flat  %
material metal  %shiny     


%%
elseif 0
    
load Sperm_Head4_Alpha0.8.mat;
fig=2;


figure(fig)
clf; set(gcf,'color','w','units','normalized');
hold on; axis equal; box off; axis off;
view(3); 

% plot sperm head.
s1= surf(xup1,xup2,xup3);
s2= surf(xdown1,xdown2,xdown3);
cc = nan(size(s1.CData,1),size(s1.CData,2),3); cc(:,:,1)=177/255; cc(:,:,2)=24/255; cc(:,:,3)=45/255; s1.CData=cc; 
cc = nan(size(s2.CData,1),size(s2.CData,2),3); cc(:,:,1)=218/255; cc(:,:,2)=207/255; cc(:,:,3)=168/255; s2.CData=cc;           

% plot sperm tail.
% Transparency parameter.
FAlpha_mat_temp = linspace(-1,0,21);
FAlpha_mat = 100.^FAlpha_mat_temp;
i_alpha_temp = 0;
for i_t = 1:size(x,2)    
    i_alpha_temp = i_alpha_temp+1;
    FAlpha = FAlpha_mat(i_alpha_temp);
    p = plot3(x(:,i_t),y(:,i_t),z(:,i_t),'color',[36 100 171]./255,'linewidth',4);    %'color',[0.8500    0.3250    0.0980]
    p.Color(4) = FAlpha;
    %PlotTailCylinder(x(:,i_t),y(:,i_t),z(:,i_t),FAlpha,fig)  
end           
shading interp
    

l1=light;
l1.Color = [1 1 1];
l1.Style = 'local';  %'infinite';  %
l1.Position = [-0.7 -0.6 -1];
l2=light;
l2.Color = [1 1 1];
l2.Style = 'local';  %'infinite';  %
l2.Position = [0.8 0 1];  
l3=light;
l3.Color = [1 1 1];
l3.Style = 'local';  %'infinite';  %
l3.Position = [1 0 1];
l4.Color = [1 1 1];
l4.Style = 'local';  %'infinite';  %
l4.Position = [1.4 0 1];
lighting gouraud   %phong  %flat  %
material metal  %shiny     
   
   

%%
else 
    
load Sperm_Head1_Alpha0.mat;
clearvars -except  xup1 xup2 xup3  xdown1 xdown2 xdown3


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
x = interp2(t0,s',Xtail,t1,s'); 
y = interp2(t0,s',Ytail,t1,s'); 
z = interp2(t0,s',Ztail,t1,s'); 
nt=size(x,2);
tRange=linspace(0,tRange(end),nt);   

load FreeSperm_Frequency.mat;  % average beat frequency: 1/HF_Freq(sp)
ncycles = tRange(end)*HF_Freq(sp);             
nclrs = round(nt/ncycles);
x = x(:,9*nclrs+1:10*nclrs)/arc_mean;
y = y(:,9*nclrs+1:10*nclrs)/arc_mean;
z = z(:,9*nclrs+1:10*nclrs)/arc_mean;
x = x-repmat(x(1,:),ns,1)+max(xup1(:));


fig=3;
figure(fig)
clf; set(gcf,'color','w','units','normalized');
hold on; axis equal; box off; axis off;
view(3); 

% plot sperm head.
s1= surf(xup1,xup2,xup3);
s2= surf(xdown1,xdown2,xdown3);
cc = nan(size(s1.CData,1),size(s1.CData,2),3); cc(:,:,1)=177/255; cc(:,:,2)=24/255; cc(:,:,3)=45/255; s1.CData=cc; 
cc = nan(size(s2.CData,1),size(s2.CData,2),3); cc(:,:,1)=218/255; cc(:,:,2)=207/255; cc(:,:,3)=168/255; s2.CData=cc;           
% plot sperm tail.
% Transparency parameter.
FAlpha_mat_temp = linspace(-1,0,nclrs);
FAlpha_mat = 100.^FAlpha_mat_temp;
i_alpha_temp = 0;  
for i_t = 1:size(x,2)    
    i_alpha_temp = i_alpha_temp+1;
    FAlpha = FAlpha_mat(i_alpha_temp);
    p = plot3(x(:,i_t),y(:,i_t),z(:,i_t),'color',[36 100 171]./255,'linewidth',4);    %'color',[0.8500    0.3250    0.0980]
    p.Color(4) = FAlpha;
    %PlotTailCylinder(x(:,i_t),y(:,i_t),z(:,i_t),FAlpha,fig)  
end           
shading interp
    


l1=light;
l1.Color = [1 1 1];
l1.Style = 'local';  %'infinite';  %
l1.Position = [-0.7 -0.6 -1];
l2=light;
l2.Color = [1 1 1];
l2.Style = 'local';  %'infinite';  %
l2.Position = [0.8 0 1];  
l3=light;
l3.Color = [1 1 1];
l3.Style = 'local';  %'infinite';  %
l3.Position = [1 0 1];
l4.Color = [1 1 1];
l4.Style = 'local';  %'infinite';  %
l4.Position = [1.4 0 1];
lighting gouraud   %phong  %flat  %
material metal  %shiny     
       
end

