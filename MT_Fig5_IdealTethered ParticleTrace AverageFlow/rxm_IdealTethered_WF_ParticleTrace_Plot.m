function rxm_IdealTethered_WF_ParticleTrace_Plot


clc

if 0 % plot sperm WF + 2 particle traces
    
% Plot sperm shape.
    load XNodes_TetheredSperm_alpha0.4_V2.mat;
    nt_WF=find(tRange==2*pi);
    Q=Nhh+Ns;   
    Xtail = Xs(Nhh+1:Q,1:nt_WF); 
    Ytail = Xs(Q+Nhh+1:2*Q,1:nt_WF); 
    Ztail = Xs(2*Q+Nhh+1:3*Q,1:nt_WF);     
    [xup1,xup2,xup3, xdown1,xdown2,xdown3] = Generate_Virtual_Head_UpDownSides;  %xup1: M1*M2.
    
    fig = 1;
    figure(fig)
    clf; set(gcf,'color','w');
    hold on; axis equal; %axis off; 
    %xlabel('x'); ylabel('y'); zlabel('z');
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////
    XL=[-0.25 1.7]; YL=[-0.3 0.5]; ZL=[-0.15 0.3]; % for tethered sperm full flow (alpha 0.4). 
    xlim([XL(1) XL(2)]); ylim([YL(1) YL(2)]); zlim([ZL(1) ZL(2)]); 
    xticks([]);yticks([]);zticks([])
    view([70 35]);
   
    [xgrid,ygrid,zgrid] = meshgrid(linspace(XL(1),XL(2),10),...
                    linspace(YL(1),YL(2),10),...
                    linspace(ZL(1),ZL(2),10));                     
CO(:,:,1) = ones(10)*0;   %ones(10)*.0; % red
CO(:,:,2) = ones(10)*0;   %ones(10)*.4470; % green
CO(:,:,3) = ones(10)*0;   %ones(10)*0.7410; % blue
surf(xgrid(:,:,1),ygrid(:,:,1),zgrid(:,:,1), CO,'FaceAlpha',0.02)% z-plane
surf(reshape(xgrid(end,:,:),[10,10]),reshape(ygrid(end,:,:),[10,10]),reshape(zgrid(end,:,:),[10,10]), CO,'FaceAlpha',0.1)%y plane
surf(reshape(xgrid(:,1,:),[10,10]),reshape(ygrid(:,1,:),[10,10]),reshape(zgrid(:,1,:),[10,10]), CO,'FaceAlpha',0.2)%x plane
shading interp


if 1    
    CO_up(:,:,1) = 0.6*ones(size(xup1)); CO_up(:,:,2) = 0*ones(size(xup1)); CO_up(:,:,3) = 0*ones(size(xup1));
    CO_down(:,:,1) = 0*ones(size(xdown1)); CO_down(:,:,2) = 0.6*ones(size(xdown1)); CO_down(:,:,3) = 0*ones(size(xdown1));
    s1= surf(xup1,xup2,xup3,CO_up);
    s2= surf(xdown1,xdown2,xdown3,CO_down);
    s1.EdgeColor='none' ;
    s2.EdgeColor='none' ;
    
    nt0 = 18;
    Xtail = [Xtail Xtail]; Ytail = [Ytail Ytail]; Ztail = [Ztail Ztail];
    clrs = colormap(winter(nt_WF));
    for i_nt = nt0:nt0+nt_WF-1
        if i_nt==nt0
            plot3(Xtail(:,i_nt),Ytail(:,i_nt),Ztail(:,i_nt),'color',clrs(i_nt-nt0+1,:),'linewidth',3) 
        else
            p1 = plot3(Xtail(:,i_nt),Ytail(:,i_nt),Ztail(:,i_nt),'color',clrs(i_nt-nt0+1,:),'linewidth',2);
            p1.Color(4) = 0.1;
        end
    end
end


%% Plot near-head particle trace (time: 100 BC).

load XNodes_TetheredSperm_alpha0.4_ParticleTrace_Head_selected_V2.mat;  %x_temp=-0.3; y_temp=0.5; z_temp=0.1;
%load XNodes_TetheredSperm_alpha0.4_ParticleTrace_Head.mat;  %x_temp=-0.3:0.1:-0.1; y_temp=0.2:0.1:0.5; z_temp=-0.1:0.1:0.1;

np = size(particle_trace,1)/3;
nt = size(particle_trace,2); 
nt_plot1 = round(nt*0.1);
nt_plot2 = round(nt*0.6);
particle_trace = particle_trace(:,nt_plot1:nt_plot2);

for i_p = 1 %1:np   
    particle_trace_temp = [particle_trace(i_p,:);  particle_trace(i_p+np,:);  particle_trace(i_p+2*np,:) ];
    plot3(particle_trace_temp(1,:),particle_trace_temp(2,:),particle_trace_temp(3,:),'k','linewidth',2)  
    plot3(particle_trace_temp(1,:),particle_trace_temp(2,:),particle_trace_temp(3,:)*0+ZL(1),'k','linewidth',1)  
end


%% Plot distal-tail particle trace (time: 100 BC).

load XNodes_TetheredSperm_alpha0.4_ParticleTrace_Tail_selected_V2.mat;  %x_temp=0.9; y_temp=0.18; z_temp=0;
%load XNodes_TetheredSperm_alpha0.4_ParticleTrace_Tail.mat; %x_temp=0.8:0.1:1.0; y_temp=0.18:0.01:0.22; z_temp=-0.1:0.1:0.1;

np = size(particle_trace,1)/3;
nt = size(particle_trace,2); 
nt_plot = round(nt*0.4);
particle_trace = particle_trace(:,1:nt_plot);

for i_p = 1:np   
    particle_trace_temp = [particle_trace(i_p,:);  particle_trace(i_p+np,:);  particle_trace(i_p+2*np,:) ];
    plot3(particle_trace_temp(1,:),particle_trace_temp(2,:),particle_trace_temp(3,:),'k','linewidth',2)  
    plot3(particle_trace_temp(1,:),particle_trace_temp(2,:),particle_trace_temp(3,:)*0+ZL(1),'k','linewidth',1)  
end


% Plot Cartesian axes.
    O = [XL(1)+0.03 YL(1)+0.03 ZL(2)-0.03]; % for PCA modes.
    Cx=[1; 0; 0].*1; 
    Cy=[0; 1; 0].*1; 
    Cz=[0; 0; -1].*1; 
    C=[Cx Cy Cz];
    Plot_Cartesian_Axes(fig,O,C)  
    camlight right
    lighting gouraud  
    
    
%%
elseif 0 % plot 1 particle trace -- particle near head.
    load XNodes_TetheredSperm_alpha0.4_ParticleTrace_Head_selected_V2.mat;  %x_temp=-0.3; y_temp=0.5; z_temp=0.1;
    %load XNodes_TetheredSperm_alpha0.4_ParticleTrace_Head.mat;  %x_temp=-0.3:0.1:-0.1; y_temp=0.2:0.1:0.5; z_temp=-0.1:0.1:0.1;

np = size(particle_trace,1)/3;
nt = size(particle_trace,2); 
nt_plot1 = round(nt*0.15);
nt_plot2 = round(nt*0.55);
particle_trace = particle_trace(:,nt_plot1:nt_plot2);

fig = 2;
figure(fig)
    clf; set(gcf,'color','w');
    hold on; axis equal; box on; %axis off; 
    %xlabel('x'); ylabel('y'); zlabel('z');
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////
    XL=[-0.16 0.45]; YL=[0.2 0.48]; ZL=[0.01 0.25]; % for tethered sperm full flow (alpha 0.4). 
    xlim([XL(1) XL(2)]); ylim([YL(1) YL(2)]); zlim([ZL(1) ZL(2)]); 
    xticks([]);yticks([]);zticks([])
    view([70 35]);
   
    [xgrid,ygrid,zgrid] = meshgrid(linspace(XL(1),XL(2),10),...
                    linspace(YL(1),YL(2),10),...
                    linspace(ZL(1),ZL(2),10));                     
CO(:,:,1) = ones(10)*0;   %ones(10)*.0; % red
CO(:,:,2) = ones(10)*0;   %ones(10)*.4470; % green
CO(:,:,3) = ones(10)*0;   %ones(10)*0.7410; % blue
surf(xgrid(:,:,1),ygrid(:,:,1),zgrid(:,:,1), CO,'FaceAlpha',0.)% z-plane
surf(reshape(xgrid(end,:,:),[10,10]),reshape(ygrid(end,:,:),[10,10]),reshape(zgrid(end,:,:),[10,10]), CO,'FaceAlpha',0.1)%y plane
surf(reshape(xgrid(:,1,:),[10,10]),reshape(ygrid(:,1,:),[10,10]),reshape(zgrid(:,1,:),[10,10]), CO,'FaceAlpha',0.2)%x plane
shading interp


for i_p = 1 %1:np   
    particle_trace_temp = [particle_trace(i_p,:);  particle_trace(i_p+np,:);  particle_trace(i_p+2*np,:) ];
    plot3(particle_trace_temp(1,:),particle_trace_temp(2,:),particle_trace_temp(3,:),'k','linewidth',2)  
    plot3(particle_trace_temp(1,:),particle_trace_temp(2,:),particle_trace_temp(3,:)*0+ZL(1),'k','linewidth',1)  
end

% Plot Cartesian axes.
    O = [XL(1)+0.03 YL(1)+0.03 ZL(2)-0.03]; % for PCA modes.
    Cx=[1; 0; 0].*0.5; 
    Cy=[0; 1; 0].*0.5; 
    Cz=[0; 0; -1].*0.5; 
    C=[Cx Cy Cz];
    Plot_Cartesian_Axes(fig,O,C)  
    camlight right
    lighting gouraud  
    
%%
elseif 1 % plot 1 particle trace -- particle near distal tail.
     
load XNodes_TetheredSperm_alpha0.4_ParticleTrace_Tail_selected_V2.mat;  %x_temp=0.9; y_temp=0.18; z_temp=0;
%load XNodes_TetheredSperm_alpha0.4_ParticleTrace_Tail.mat; %x_temp=0.8:0.1:1.0; y_temp=0.18:0.01:0.22; z_temp=-0.1:0.1:0.1;

np = size(particle_trace,1)/3;
nt = size(particle_trace,2); 
nt_plot1 = round(nt*0.02);
nt_plot2 = round(nt*0.42);
particle_trace = particle_trace(:,nt_plot1:nt_plot2);

fig=3;
figure(fig)
    clf; set(gcf,'color','w');
    hold on; axis equal; box on; %axis off; 
    %xlabel('x'); ylabel('y'); zlabel('z');
    %////////////////////////////////////////////////////////////////////////////////////////////////////////////
    XL=[0.95 1.7]; YL=[-0. 0.3]; ZL=[-0.03 0.3]; % for tethered sperm full flow (alpha 0.4). 
    xlim([XL(1) XL(2)]); ylim([YL(1) YL(2)]); zlim([ZL(1) ZL(2)]); 
    xticks([]);yticks([]);zticks([])
    view([70 35]);
    
    
    [xgrid,ygrid,zgrid] = meshgrid(linspace(XL(1),XL(2),10),...
                    linspace(YL(1),YL(2),10),...
                    linspace(ZL(1),ZL(2),10));                     
CO(:,:,1) = ones(10)*0;   %ones(10)*.0; % red
CO(:,:,2) = ones(10)*0;   %ones(10)*.4470; % green
CO(:,:,3) = ones(10)*0;   %ones(10)*0.7410; % blue
surf(xgrid(:,:,1),ygrid(:,:,1),zgrid(:,:,1), CO,'FaceAlpha',0.)% z-plane
surf(reshape(xgrid(end,:,:),[10,10]),reshape(ygrid(end,:,:),[10,10]),reshape(zgrid(end,:,:),[10,10]), CO,'FaceAlpha',0.1)%y plane
surf(reshape(xgrid(:,1,:),[10,10]),reshape(ygrid(:,1,:),[10,10]),reshape(zgrid(:,1,:),[10,10]), CO,'FaceAlpha',0.2)%x plane
shading interp
    
for i_p = 1:np   
    particle_trace_temp = [particle_trace(i_p,:);  particle_trace(i_p+np,:);  particle_trace(i_p+2*np,:) ];
    plot3(particle_trace_temp(1,:),particle_trace_temp(2,:),particle_trace_temp(3,:),'k','linewidth',2)  
    plot3(particle_trace_temp(1,:),particle_trace_temp(2,:),particle_trace_temp(3,:)*0+ZL(1),'k','linewidth',1)  
end


% Plot Cartesian axes.
    O = [XL(1)+0.03 YL(1)+0.03 ZL(2)-0.03]; % for PCA modes.
    Cx=[1; 0; 0].*0.6; 
    Cy=[0; 1; 0].*0.6; 
    Cz=[0; 0; -1].*0.6; 
    C=[Cx Cy Cz];
    Plot_Cartesian_Axes(fig,O,C)  
    camlight right
    lighting gouraud  
end
