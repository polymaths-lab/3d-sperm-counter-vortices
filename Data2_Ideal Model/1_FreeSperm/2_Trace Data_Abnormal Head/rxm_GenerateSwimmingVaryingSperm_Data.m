function rxm_GenerateSwimmingVaryingSperm_Data

clc


%% Varying setting-up for various virtual sperm cases.
% Various head sizes.
a1_mat = [2.0/45  10/45  20/45  30/45]; 
a2_mat = [1.6/45  8/45   16/45  24/45];
a3_mat = [1.0/45  5/45   10/45  15/45];
% Various waveform planarity (2D/ 3D).
alpha_mat = [0 0.2 0.4 0.6 0.8];
% Various boundary conditions.
domain_mat = {'i'; 'h'};   %domain='i'--> Stokeslet; domain='h'--> Blakelet.



%% Consistent setting-up.
nbeats=1;%31; 
dt=2*pi*0.05;
tRange=0:dt:2*pi*nbeats;
nns=30; 
s = linspace(0,1,nns);
[S,T]=ndgrid(s,tRange);

args.phase=0;
args.k=2*pi;

B{1}=RotationMatrix(0*pi/3,3);
b10{1}=B{1}(:,1);
b20{1}=B{1}(:,2);


bound_dis = 0.2-1.0/45; %uniformed boundary distance: from the boundary to the initial moment of the sperm head surface's bottom point.
boundary=[];

% discretisation parameters - number of points
swimmer{1}.model.ns=40;
swimmer{1}.model.nh=4;
swimmer{1}.model.Ns=100;
swimmer{1}.model.Nh=10;

epsilon=0.25/45;
blockSize=0.2;

HC_dis = nan(length(a1_mat), length(alpha_mat), length(domain_mat)); % sperm head center swimming distance.

%% Generate swimmer and solve swimming problem.

for i_head = 1:length(a1_mat)
    for i_alpha = 1:length(alpha_mat)
        for i_bound = 1:length(domain_mat)
            % Prescribe specific setting-up for each virtual sperm case.
            xyWaveFn=@DKActSpermWave_modified;
            x00{1}=[0;0;bound_dis+a3_mat(i_head)];
            swimmer{1}.model.F=ConstructInterpolantFromxyForm_modified(S,T,xyWaveFn,alpha_mat(i_alpha),args);
            swimmer{1}.fn=@SpermModelGI;
            swimmer{1}.model.a1=a1_mat(i_head);
            swimmer{1}.model.a2=a2_mat(i_head);
            swimmer{1}.model.a3=a3_mat(i_head);
            domain=domain_mat{i_bound}; 

            
            if 1 % Plot the waveform to check if the prescribed sperm model is correct.
                clrs = summer(length(tRange));
                figure(1)
                clf; hold on;set(gcf,'color','w'); view(3); box on; axis equal;
                for i_t=1:length(tRange)
                    xlabel('$x$','interpreter','latex'); ylabel('$y$','interpreter','latex'); zlabel('$z$','interpreter','latex');
                    [xis,~,~]  = swimmer{1}.fn(tRange(i_t),swimmer{1}.model);
                    [x1,x2,x3] = ExtractComponents(xis);
                    plot3(x1,x2,x3,'.','MarkerSize',20,'color',clrs(i_t,:));
                end
            end

            

            tic
            fprintf('starting solver\n')
            [t,z]=SolveSwimmingTrajectoryAndForces(x00{1},b10{1},b20{1},tRange,swimmer{1},...
            boundary,epsilon,domain,blockSize);
            solveTime = toc;
            fprintf('CPU time taken = %f\n',solveTime)
            
            
            HC_dis_temp = arclength(z(:,1),z(:,2),z(:,3)); 
            HC_dis(i_head,i_alpha,i_bound) = HC_dis_temp(end); 
        end
    end
end


%% Data.

% relative deviation of the sperm head center swimming distances between with-boundary and no-boundary cases.
Dis_RelDev = (  HC_dis(:,:,2)-HC_dis(:,:,1) )./ HC_dis(:,:,1) 


Head_a1 = a1_mat;
Head_a2 = a2_mat;
Head_a3 = a3_mat;

WF_alpha = alpha_mat;


save('DisDev_VarySperm_31BC.mat','Head_a1','Head_a2','Head_a3','WF_alpha','Dis_RelDev');

