function rxm_IdealWF_data

clc

nbeats=1; 
tRange=[0 2*pi*nbeats];
nns=30; 

dt=2*pi*0.05;  %dt=2*pi*0.01;%

% Generate swimmer ----------------------------------
% waveform
xyWaveFn=@DKActSpermWave;
args.phase=0;
args.k=2*pi;  
s = linspace(0,1,nns);
t = tRange(1):dt:tRange(2);
[S,T]=ndgrid(s,t);
x00{1}=[0;0;0.2];
B{1}=RotationMatrix(0*pi/3,3);
b10{1}=B{1}(:,1);
b20{1}=B{1}(:,2);
swimmer{1}.model.F=ConstructInterpolantFromxyForm(S,T,xyWaveFn,args);
swimmer{1}.fn=@SpermModelGI;

% discretisation parameters - number of points
swimmer{1}.model.ns=40;
swimmer{1}.model.nh=4;
swimmer{1}.model.Ns=100;
swimmer{1}.model.Nh=10;
% semi-axes
swimmer{1}.model.a1=2.0/45;
swimmer{1}.model.a2=1.6/45;
swimmer{1}.model.a3=1.0/45;

% Panel a ----------------------------------
figure(1);clf;hold on;
nth=10;
nphi=20;
a=1;
[x3,x2,x1]=GenerateSphereSurfaceForVisualisation(nth,nphi,a);
x1=x1*swimmer{1}.model.a1;
x2=x2*swimmer{1}.model.a2;
x3=x3*swimmer{1}.model.a3;

surf(x1,x2,x3,0*x3);shading flat;light;

clrs = summer(length(t));
for nt=1:length(t)
    [~,~,xis]  = swimmer{1}.fn(t(nt),swimmer{1}.model);
    [x1,x2,x3] = ExtractComponents(xis);
    plot3(x1,x2,x3,'.','color',clrs(nt,:));
    x(:,nt) = x1(601:end)-swimmer{1}.model.a1;
    y(:,nt) = x2(601:end);
    z(:,nt) = x3(601:end);
end
axis equal;axis on;

%% Save data


save('IdealWF_alpha0.4.mat','x','y','z','t','s');
%save('IdealWF_alpha0.6_Tfine.mat','x','y','z','t','s');

