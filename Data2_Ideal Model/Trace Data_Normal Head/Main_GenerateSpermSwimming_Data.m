function Main_GenerateSpermSwimming_Data

clc


%% sperm waveform (image)


nbeats=1;
tRange=[0 2*pi*nbeats];
nns=30; 
dt=2*pi*0.05;

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



figure(1)
clf; hold on;
axis on;axis equal;
clrs = summer(length(t));
for nt=1:length(t)
    [xis,~,~]  = swimmer{1}.fn(t(nt),swimmer{1}.model);
    [x1,x2,x3] = ExtractComponents(xis);
    plot3(x1,x2,x3,'.','color',clrs(nt,:));
end





%% Calculate sperm swimming.


% Define waveform.
nns=30; 
s = linspace(0,1,nns);
nbeats=20;
tRange_lim=[0 2*pi*nbeats];
dt=2*pi*0.05;
tRange= tRange_lim(1):dt:tRange_lim(2);
[S,T]=ndgrid(s,tRange);
xyWaveFn=@DKActSpermWave;
args.phase=0;
args.k=2*pi;
swimmer{1}.model.F=ConstructInterpolantFromxyForm(S,T,xyWaveFn,args);

% discretisation parameters - number of points
swimmer{1}.model.ns=40;
swimmer{1}.model.nh=4;
swimmer{1}.model.Ns=100;
swimmer{1}.model.Nh=10;
% semi-axes
swimmer{1}.model.a1=2.0/45;
swimmer{1}.model.a2=1.6/45;
swimmer{1}.model.a3=1.0/45;

swimmer{1}.fn=@SpermModelGI;
x00{1}=[0;0;0.2];
B{1}=RotationMatrix(0*pi/3,3);
b10{1}=B{1}(:,1);
b20{1}=B{1}(:,2);

% boundary
boundary=[];

% numerical parameters
epsilon=0.25/45;
domain='i';
blockSize=0.2;


tic
fprintf('starting solver\n')
[t,z]=SolveSwimmingTrajectoryAndForces(x00{1},b10{1},b20{1},tRange,swimmer{1},...
   boundary,epsilon,domain,blockSize);
solveTime = toc;

fprintf('CPU time taken = %f\n',solveTime)



%% Figure to check the lab trajectory. 

figure(2);clf;hold on;view(3)
Nsw=length(swimmer);
x0 = cell(1,Nsw);
for n=1:Nsw
    x0{n}=z(:,n:Nsw:n+2*Nsw);
    plot3(x0{n}(:,1),x0{n}(:,2),x0{n}(:,3),'k');
end
hx=xlabel('\(x_1\) (flagellar lengths)','interpreter','latex');
hy=ylabel('\(x_2\) (flagellar lengths)','interpreter','latex');
hz=zlabel('\(x_3\) (flagellar lengths)','interpreter','latex');
%axis equal;
%xlim([-1.5,1.5]);
%ylim([-1.5,1.5]);
box on;
set(gca,'tickdir','out');


%% Save data.

head_x0=z(:,1:3); 
head_tangent=z(:,4:6);
head_normal=z(:,7:9);
H = z(:,10:end);

ns=swimmer{1}.model.ns;
Ns=swimmer{1}.model.Ns;
nt = length(t);

xb = [];
Xb = [];
nb = 0;
Nb = 0;


for i_nt = 1:nt
    i_nt
    [xi,~,Xi]=swimmer{1}.fn(t(i_nt),swimmer{1}.model);
    nhh = length(xi)/3-ns;
    Nhh = length(Xi)/3-Ns;
  
    x0=z(i_nt,1:3);
    b1=z(i_nt,4:6);
    b2=z(i_nt,7:9);
    b3=cross(b1,b2);
    B=[b1(:) b2(:) b3(:)];
    
    % Flagellar force and quadrature points.
    xx0=ApplyRotationMatrix(B,xi); 
    xs(:,i_nt)=TranslatePoints(xx0,x0);
    Xx0=ApplyRotationMatrix(B,Xi);
    Xs(:,i_nt)=TranslatePoints(Xx0,x0);

    % Get head index.
    if i_nt==1
        f=nan(3*(nhh+ns+nb),nt);
       
        [xsperm1_temp,~,xsperm3_temp] = ExtractComponents(xi);
        xh1 = xsperm1_temp(1:nhh);
        xh3 = xsperm3_temp(1:nhh);
        Ind_head_up = find(xh3>=0);  
        Ind_head_down = find(xh3<0);
        [~,Ind_head_neck] = max(xh1);
        
        [Xsperm1_temp,~,Xsperm3_temp] = ExtractComponents(Xi);
        Xh1 = Xsperm1_temp(1:Nhh);
        Xh3 = Xsperm3_temp(1:Nhh);
        Ind_Head_up = find(Xh3>=0);  
        Ind_Head_down = find(Xh3<0);
        [~,Ind_Head_neck] = max(Xh1);
    end
    
    % force distribution
    f_temp=ExtractForceDistributionFromTimeIntegral(tRange,H,tRange(i_nt)); 
    f(:,i_nt)=f_temp(:);
   
end

save(['XNodes_Lab_alpha0.6.mat'],'xs','Xs','xb','Xb','f','nhh','Nhh','ns','Ns','nb','Nb',...
     'Ind_head_up','Ind_head_down','Ind_head_neck','Ind_Head_up','Ind_Head_down','Ind_Head_neck',...
     'head_x0','head_tangent','head_normal','tRange');
