function rxm_NumRecBlakelet_MobilityXNodes_Data

% Under the numerical framework of regularized Stokeslet method, this 
% function numerically reconstructs the experimentally observed sperm 
% mobility through PCA-reconstructed experimental waveform. 
% Both with-one-boundary and without boundary cases are considered in our
% simulations, but this function is merely for the no-boundary case.


clc

% near-boundary sperm.
sp=8; 

load 03_Head_Fixed_Frame_Method_Rotations_Traces_raw_5_Microns_Difraction_head_spin_rxmHeadNBT.mat
Xdata = X{sp};
Ydata = Y{sp};
Zdata = Z{sp};
ns=size(Xdata,1);
nt=size(Xdata,2);
 
load FreeSperm_Frequency.mat; % read the beat frequency (HF WF frequency) to encode the WF color
Freq_WF = HF_Freq(sp);
% (time in seconds)  
dt = 1/90; % Each frame is taken every (1/90) seconds => 1.00 s corresponds to 90 video data frames (framerate was 90fps).
tRange = dt*(nt-1);
ncycle=Freq_WF*tRange;   
   


%% Step1: reconstruct the observed sperm waveform using PCA.

%======================================================================================================================================================
% Extract PCA components.
    [~,C,B]=Extract_XYZ_PCA(Xdata,Ydata,Zdata); %output includes 1-22 PCA modes, nmode=22.
    % Output:
    % D: sorted eigenvalues of the decomposed PCA components.
    % C: sorted eigenvectors of the decomposed PCA components. ns*((1+nmode)*3).
    % B: time coefficients of the sorted eigenvectors. nt*nmode.
    nmode = size(B,2);
    for i_nt=1:size(B,1)        
        Xrec122(:,:,i_nt)= C(:,1:3); 
        for i_nV =1:nmode
            Xrec122(:,:,i_nt)= Xrec122(:,:,i_nt)+ B(i_nt,i_nV)*C(:, 1+i_nV*3 : 3+i_nV*3);% superimpose PCA modes 0th-22nd (1st-22nd)   
        end      
        xrec122(:,i_nt)=Xrec122(:,1,i_nt);
        yrec122(:,i_nt)=Xrec122(:,2,i_nt);
        zrec122(:,i_nt)=Xrec122(:,3,i_nt);   
    end
 
%save(['WF_PCA_sp',num2str(sp),'.mat'],'D','C','B','T');




%======================================================================================================================================================
% Plot the reconstructed waveform to check the PCA implementation.
% Font options for printing figures
fs=30;
fn='times';

nclrs = round(nt/ncycle);
clrs = colormap(parula(nclrs));

figure(1)
clf;
xmax=70;ymax=50;zmin=-40;
%raw waveform
subplot(1,2,1)
hold on;view(3);axis equal;
for i_nt=1:nt
    plot3(Xdata(:,i_nt),Ydata(:,i_nt),Zdata(:,i_nt),'color',clrs(mod(i_nt,nclrs)+1,:))    
    % projections
    plot3(xmax*ones(ns,1),Ydata(:,i_nt),Zdata(:,i_nt),'color',[0.8,0.8,0.8],'linewidth',0.1)
    plot3(Xdata(:,i_nt),ymax*ones(ns,1),Zdata(:,i_nt),'color',[0.8,0.8,0.8],'linewidth',0.1)
    plot3(Xdata(:,i_nt),Ydata(:,i_nt),zmin*ones(ns,1),'color',[0.8,0.8,0.8],'linewidth',0.1)   
end
xlim([-5 xmax]); ylim([-30 ymax]);zlim([zmin 25]);
set(gca,'fontsize',fs); set(gca,'fontname',fn);
hx=xlabel('x (\mum)');hy=ylabel('y (\mum)');hz=zlabel('z (\mum)');
set(hx,'fontsize',fs);  set(hx,'fontname',fn);
set(hy,'fontsize',fs);  set(hy,'fontname',fn);
set(hz,'fontsize',fs);  set(hz,'fontname',fn);
title('Raw','fontsize',fs)
%
% reconstructed WF using 1-22 PCA modes
subplot(1,2,2)
hold on;view(3);axis equal;
for i_nt=1:nt   
    plot3(xrec122(:,i_nt),yrec122(:,i_nt),zrec122(:,i_nt),'color',clrs(mod(i_nt,nclrs)+1,:))    
    %projection
    plot3(xmax*ones(ns,1),yrec122(:,i_nt),zrec122(:,i_nt),'color',[0.8,0.8,0.8],'linewidth',0.1)
    plot3(xrec122(:,i_nt),ymax*ones(ns,1),zrec122(:,i_nt),'color',[0.8,0.8,0.8],'linewidth',0.1)
    plot3(xrec122(:,i_nt),yrec122(:,i_nt),zmin*ones(ns,1),'color',[0.8,0.8,0.8],'linewidth',0.1)   
end
xlim([-5 xmax]); ylim([-30 ymax]);zlim([zmin 25]);
set(gca,'fontsize',fs); set(gca,'fontname',fn);
hx=xlabel('x (\mum)');hy=ylabel('y (\mum)');hz=zlabel('z (\mum)');
set(hx,'fontsize',fs);  set(hx,'fontname',fn);
set(hy,'fontsize',fs);  set(hy,'fontname',fn);
set(hz,'fontsize',fs);  set(hz,'fontname',fn);
title('WF using PCA modes 1-22','fontsize',fs)%title('2 PCA modes','fontsize',fs)


%spt= suptitle(['Waveform in the head-fixed frame: sp#',num2str(sp)]);
%set(spt,'fontsize',fs+2)
%print(1,'-dpng','-r600',['hfWF_Reconstruct_via_PCA_sp',num2str(sp),'.png']);




%======================================================================================================================================================
% Spline-smooth PCA coefficients 'B(t)', via function 'csaps'.

NT=1:1:nt; % timepoint sequence
t_pca = linspace(0,tRange,100*round(ncycle));
NT_pca = linspace(1,nt,100*round(ncycle));% timepoint sequence

pp1=0.99; %'pp=1' means interpolation without cubic spline.
for i_nV=1:size(B,2)
    B_ss(:,i_nV)=fnval(csaps(NT,B(:,i_nV),pp1),NT_pca);
end
%save(['WF_PCA_sp',num2str(sp),'.mat'],'B_splinesmooth_BCs','w_fft','t_pca','-append');



%======================================================================================================================================================
% Plot the cubic-spline-smoothed coefficients.
B1=B(:,1);
Bss1 = B_ss(:,1);

figure(2)
clf;
subplot(211)
hold on;
scatter(NT,B1,15,'r','filled')
plot(NT_pca,Bss1,'b')
set(gca,'fontsize',fs); set(gca,'fontname',fn);
hy=ylabel('PCA coefficient');
set(hy,'fontsize',fs);  set(hy,'fontname',fn);
title('Over a longer time','fontsize',fs)

subplot(212)
hold on;
scatter(NT(1:nclrs),B1(1:nclrs),15,'r','filled')
plot(NT_pca(1:100),Bss1(1:100),'b')
set(gca,'fontsize',fs); set(gca,'fontname',fn);
hx=xlabel('time point sequence');
hy=ylabel('PCA coefficient');
set(hx,'fontsize',fs);  set(hx,'fontname',fn);
set(hy,'fontsize',fs);  set(hy,'fontname',fn);
legend({'B1(t)','B1, cubic spline smooth: pp 0.99'},'location','best');
title(['Over 1-beat-cycle time'],'fontsize',fs)

%spt = suptitle(['Cubic spline smoothing of PCA coefficient B1, WF of sp#',num2str(sp)]);
%set(spt,'fontsize',fs+2);




%======================================================================================================================================================
% Generate the PCA-reconstructed WF, using the smoothed coefficients B(t).

ntt=size(B_ss,1);
for i_nt=1:ntt
    Xrec122_Bss(:,:,i_nt)= C(:,1:3); 
    for i_nV =1:size(B_ss,2)
        Xrec122_Bss(:,:,i_nt)= Xrec122_Bss(:,:,i_nt)+ B_ss(i_nt,i_nV)*C(:, 1+i_nV*3 : 3+i_nV*3);% superimpose PCA modes 0th-22nd (1st-22nd)
    end
    xrec122_Bss(:,i_nt)=Xrec122_Bss(:,1,i_nt);
    yrec122_Bss(:,i_nt)=Xrec122_Bss(:,2,i_nt);
    zrec122_Bss(:,i_nt)=Xrec122_Bss(:,3,i_nt);    
end

%save(['WF_PCA_sp',num2str(sp),'.mat'],'D','C','B','T','B_splinesmooth_BCs','w_fft','t_pca','xm122_Bss','ym122_Bss','zm122_Bss');




%======================================================================================================================================================
% Normalize the PCA-reconstructed WF.

clear X Y Z

X=xrec122_Bss;
Y=yrec122_Bss;
Z=zrec122_Bss;
ns=size(X,1);nt=size(X,2);
a=arclength(X,Y,Z);

X_n=nan(ns,nt);Y_n=X_n;Z_n=X_n; % normalized coordinates
for i_nt=1:nt
    X_n(:,i_nt)=X(:,i_nt)./a(end,i_nt);
    Y_n(:,i_nt)=Y(:,i_nt)./a(end,i_nt);
    Z_n(:,i_nt)=Z(:,i_nt)./a(end,i_nt);
end
X_n=X_n-repmat(X_n(1,:),ns,1);
Y_n=Y_n-repmat(Y_n(1,:),ns,1);
Z_n=Z_n-repmat(Z_n(1,:),ns,1);


S=linspace(0,1,ns); % normalized arclength
s=linspace(0,1,40);
%The following interpolation is for the coarse grid of sperm tail coordinates at the step of 'nearest neighbour matrix'.
x_n = interp2(t_pca,S',X_n,t_pca,s'); 
y_n = interp2(t_pca,S',Y_n,t_pca,s'); 
z_n = interp2(t_pca,S',Z_n,t_pca,s'); 

X_pca122=X_n;  Y_pca122=Y_n;  Z_pca122=Z_n;
x_pca122=x_n;  y_pca122=y_n;  z_pca122=z_n;

%save(['NormalizedWF_PCA_Bsplinesmooth_sp',num2str(sp),'.mat'],'X_pca122','Y_pca122','Z_pca122','x_pca122','y_pca122','z_pca122','t_pca','w_fft');%'-append')%,



%======================================================================================================================================================
% Plot the smoothed PCA-reconstructed waveform and its normalized counterpart for checking purpose.

figure(3)
set(gcf,'color','w','units','normalized','outerposition',[0 0 1 1]) % make it white background
clf;

ncycle = t_pca(end)*Freq_WF;
nclrs = round(nt/ncycle);
clrs=parula(nclrs);

subplot(121)
hold on;
view(3); axis equal;
xmax=60;ymax=40;zmin=-30;
for i=1:nclrs 
    plot3(X(:,i),Y(:,i),Z(:,i),'color',clrs(mod(i,nclrs)+1,:))
    % projection
    plot3(xmax*ones(size(X(:,i))),Y(:,i),Z(:,i),'color',[0.8,0.8,0.8],'linewidth',0.1)
    plot3(X(:,i),ymax*ones(size(X(:,i))),Z(:,i),'color',[0.8,0.8,0.8],'linewidth',0.1)
    plot3(X(:,i),Y(:,i),zmin*ones(size(X(:,i))),'color',[0.8,0.8,0.8],'linewidth',0.1)
end
set(gca,'fontsize',fs); set(gca,'fontname',fn);
hx=xlabel('x (\mum)');hy=ylabel('y (\mum)');hz=zlabel('z (\mum)');
set(hx,'fontsize',fs);  set(hx,'fontname',fn);
set(hy,'fontsize',fs);  set(hy,'fontname',fn);
set(hz,'fontsize',fs);  set(hz,'fontname',fn);
title(['original'],'fontsize',fs)


subplot(122)
hold on;
view(3); axis equal;
zlim([-0.7 0.4]);
xmax=1.5;ymax=1;zmin=-0.7;
for i=1:nclrs 
    plot3(X_n(:,i),Y_n(:,i),Z_n(:,i),'color',clrs(mod(i,nclrs)+1,:))
    % projection
    plot3(xmax*ones(size(X_n(:,i))),Y_n(:,i),Z_n(:,i),'color',[0.8,0.8,0.8],'linewidth',0.1)
    plot3(X_n(:,i),ymax*ones(size(X_n(:,i))),Z_n(:,i),'color',[0.8,0.8,0.8],'linewidth',0.1)
    plot3(X_n(:,i),Y_n(:,i),zmin*ones(size(X_n(:,i))),'color',[0.8,0.8,0.8],'linewidth',0.1)
end
set(gca,'fontsize',fs); set(gca,'fontname',fn);
hx=xlabel('x');hy=ylabel('y');hz=zlabel('z');
set(hx,'fontsize',fs);  set(hx,'fontname',fn);
set(hy,'fontsize',fs);  set(hy,'fontname',fn);
set(hz,'fontsize',fs);  set(hz,'fontname',fn);
title(['normalised'],'fontsize',fs)

%suptitle(['WF in HOV frame, sp#',num2str(sp)])

%print(3,'-dpng','-r600',['NormalizedWF_HF_MyPCAlimitcycle_sp',num2str(sp),'.png']);




%% Step2: calculate the lab trajectory using the PCA-reconstructed waveform.

%======================================================================================================================================================
% Reconstruction calculation of the observed sperm motility.

% Head position.==============================================
% discretisation parameters - number of points
model.ns=40;
model.nh=4;
model.Ns=100;
model.Nh=10;
% semi-axes
model.a1=2.0/45;
model.a2=1.6/45;
model.a3=1.0/45;

xh=GenerateSphereDiscr(model.nh,1);
[xh1,xh2,xh3]=ExtractComponents(xh);
xh1=model.a1*xh1;
xh2=model.a2*xh2;
xh3=model.a3*xh3;
xhh = [xh1;xh2;xh3];
Ind_head_up = find(xh3>=0);  
Ind_head_down = find(xh3<0);
[~,Ind_head_neck] = max(xh1);
% head is stationary in body frame
vh1=0*xh1;
vh2=0*xh2;
vh3=0*xh3;
vh = [vh1;vh2;vh3];
% generate head position
Xh=GenerateSphereDiscr(model.Nh,1);
% size(Xh)=[1800,1]
[Xh1,Xh2,Xh3]=ExtractComponents(Xh);
Xh1=model.a1*Xh1;% size(Xh1)=[600,1]
Xh2=model.a2*Xh2;
Xh3=model.a3*Xh3;
Xhh = [Xh1;Xh2;Xh3];% size(Xhh)=[1800,1]
Ind_Head_up = find(Xh3>=0);  
Ind_Head_down = find(Xh3<0);
[~,Ind_Head_neck] = max(Xh1);


% Tail position.===================================================
%load NormalizedWF_PCA_Bsplinesmooth_sp23.mat;
Xt=X_pca122;Yt=Y_pca122;Zt=Z_pca122;
xt=x_pca122;yt=y_pca122;zt=z_pca122;
% Ns*nt arrays, 'Ns' is the discretisation parameter for tail, Ns=100.
% 'Ns' for fine grid - position only
Xtt = [Xt+ones(size(Xt))*model.a1;Yt;Zt]; % (3*Ns)*nt array
% ns*nt arrays, 'ns' is the discretisation parameter for tail, ns=40.
% 'ns' for coarse grid - position and velocity
xtt=[xt+ones(size(xt))*model.a1;yt;zt]; % (3*ns)*nt array


% Initial conditions for the swimming sperm.===================================================
% In the function 'RegBlakelet', the boundary position is set as Z=0 as
% default, so the initial Z position of sperm head center represents the 
% height of boundary.
bound_h=0.5;  %0.2;
x00=[0;0;bound_h];
load 03_Head_Fixed_Frame_Method_Rotations_Traces_raw_5_Microns_Difraction_head_spin_rxmHeadNBT.mat;
b10=head_tangent{sp}(:,1);
b20=head_normal{sp}(:,1);  
clear X Y Z


% Boundary setting.===================================================
boundary=[];    
   
% Numerical parameters.===================================================
epsilon=0.25/45;
%domain='i';  %Stokeslet
domain='h';    %Blakelet
blockSize=0.2;
% Non-dimensionalize the time.
ncycle = t_pca(end)*Freq_WF;
tRange_nd = ncycle*2*pi; %non-dimensional tRange.
t_nd = t_pca./t_pca(end)*tRange_nd;  %non-dimensional time sequence.

tic
fprintf('starting solver\n')
[t,z]=SolveSwimmingTrajectoryAndForces_WF(x00,b10,b20,t_nd,t_nd,Xhh,xhh,vh,Xtt,xtt,...
   boundary,epsilon,domain,blockSize);
% Note that t_nd=t'.
solveTime = toc;
fprintf('CPU time taken = %f\n',solveTime)

%save(['LabTraj_PCAm122_sp',num2str(sp),'_With1Bound_H',num2str(bound_h),'L.mat'],'t','z');





%======================================================================================================================================================
% Plot the lab trajectory.

figure(4)
clf;hold on;
plot3(z(:,1),z(:,2),z(:,3))

daspect([1 1 1])
view(3)
xlabel('x')
ylabel('y')
zlabel('z')


title(['Lab trajectory of sp#',num2str(sp),' (PCA-reconstructed WF: B(t) Spline smooth, mode 1-22)'])
%print(4,'-dpng','-r600','LabTraj_WF1.png');






%% Step3: save XNodes data for quadrature and force points.

head_x0 = z(:,1:3);
head_tangent = z(:,4:6);
head_normal = z(:,7:9);
H = z(:,10:end); % force integral
nt = length(t_nd); %Because t_nd=t'.
clear z t B

for i_nt=1:nt
    b1 = head_tangent(i_nt,:); 
    b2 = head_normal(i_nt,:); 
    b3=cross(b1,b2);
    B=[b1(:) b2(:) b3(:)];
    x0 = head_x0(i_nt,:)';
    t0 = t_nd(i_nt); 
          
    if i_nt==1
        nhh = length(xhh)/3;
        Nhh = length(Xhh)/3;
        ns = size(xtt,1)/3;
        Ns = size(Xtt,1)/3;
        xs=nan(3*(nhh+ns),nt);
        Xs=nan(3*(Nhh+Ns),nt);
        f=nan(3*(nhh+ns),nt);
    end
    
    % Input: head and tail nodes from the normalized WF, where head center
    % is at [0 0 0].
    [xs_temp,Xs_temp] = GetQuadForceGrid(B,x0,Xhh,xhh,Xtt,xtt,t_nd,t0);
    xs(:,i_nt) = xs_temp;
    Xs(:,i_nt) = Xs_temp;
    
    % force distribution
    f_temp=ExtractForceDistributionFromTimeIntegral(t_nd,H,t0); 
    f(:,i_nt)=f_temp(:); 
end



save(['XNodes_sp',num2str(sp),'_BoundBlakelet_H',num2str(bound_h),'.mat'],'xs','Xs','f','nhh','Nhh','ns','Ns',...
     'Ind_head_up','Ind_head_down','Ind_head_neck','Ind_Head_up','Ind_Head_down','Ind_Head_neck',...
     'head_x0','head_tangent','head_normal','t_nd','Freq_WF');

 