function rxm_GetSpermShape_data

clc

nbeats=1;% 
tRange=[0 2*pi*nbeats];
nns=100; 

dt=2*pi*0.05;

% Generate swimmer ----------------------------------
% waveform
xyWaveFn=@DKActSpermWave;
args.phase=0;
args.k=2*pi;  
s = linspace(0,1,nns);
t = tRange(1):dt:tRange(2);
[S,T]=ndgrid(s,t);
%x00{1}=[0;0;0.2];
%B{1}=RotationMatrix(0*pi/3,3);
%b10{1}=B{1}(:,1);
%b20{1}=B{1}(:,2);
swimmer{1}.model.F=ConstructInterpolantFromxyForm(S,T,xyWaveFn,args);
swimmer{1}.fn=@SpermModelGI;

% discretisation parameters - number of points
swimmer{1}.model.ns=40;
swimmer{1}.model.nh=4;
swimmer{1}.model.Ns=100;
swimmer{1}.model.Nh=10;
% semi-axes
if 0
    swimmer{1}.model.a1=2.0/45;
    swimmer{1}.model.a2=1.6/45;
    swimmer{1}.model.a3=1.0/45;
else
   swimmer{1}.model.a1=30/45;
   swimmer{1}.model.a2=24/45;
   swimmer{1}.model.a3=15/45; 
end

figure(2) 
clf;
hold on; view(3)

if 0  %head plot.
nth=10;
nphi=20;
a=1;
[x3,x2,x1]=GenerateSphereSurfaceForVisualisation(nth,nphi,a);
x1=x1*swimmer{1}.model.a1;
x2=x2*swimmer{1}.model.a2;
x3=x3*swimmer{1}.model.a3;
surf(x1,x2,x3,0*x3);shading flat;light;
end


% Generate sperm head.
    %head:up+down£¬different colors,surf
    nth=20;%10;
    nphi=40;%20;
    a=1;
    M=nth*nphi;
    %head: up side
    [xup1,xup2,xup3]=GenerateSphereSurfaceForVisualisation_up(nth,nphi,a);
    [M1 M2]=size(xup1);
    xup1=xup1*swimmer{1}.model.a1;
    xup2=xup2*swimmer{1}.model.a2;
    xup3=xup3*swimmer{1}.model.a3;
    %head: down side
    [xdown1,xdown2,xdown3]=GenerateSphereSurfaceForVisualisation_down(nth,nphi,a);
    xdown1=xdown1*swimmer{1}.model.a1;
    xdown2=xdown2*swimmer{1}.model.a2;
    xdown3=xdown3*swimmer{1}.model.a3;
    CO_up(:,:,1) = 0.6*ones(size(xup1)); CO_up(:,:,2) = 0*ones(size(xup1)); CO_up(:,:,3) = 0*ones(size(xup1));
    CO_down(:,:,1) = 0*ones(size(xdown1)); CO_down(:,:,2) = 0.6*ones(size(xdown1)); CO_down(:,:,3) = 0*ones(size(xdown1));
    if 0
    s1= surf(xup1,xup2,xup3,CO_up);
    s2= surf(xdown1,xdown2,xdown3,CO_down);%,source_ligth,k
    s1.EdgeColor='none' ;
    s2.EdgeColor='none' ;
    end
    
    
    
clrs = summer(length(t));
for nt=1:length(t)
    [~,~,xis]  = swimmer{1}.fn(t(nt),swimmer{1}.model);
    [x1,x2,x3] = ExtractComponents(xis);
    x(:,nt) = x1(601:end);
    y(:,nt) = x2(601:end);
    z(:,nt) = x3(601:end);
    plot3(x1,x2,x3,'color',clrs(nt,:),'linewidth',3);   
    %plot3(x,y,z,'k','linewidth',3);    
end
axis equal;axis on;

%return
%% Save data

save('Sperm_Head4_Alpha0.8.mat','x','y','z','t','s','xup1','xup2','xup3','xdown1','xdown2','xdown3');

