function rxm_IdealTetheredSperm_FlowParticleTrace_Data
% This function calculates the particle trace in the flow field induced by the idealized sperm, where the sperm head center is located at the origin point.

clc

%p1_x0 = [-0.2; 0.5; 0.1]; %离精子太远了，和参考文献不符。
%p2_x0 = [0.8; 0.2; 0]; %离精子稍远，与参考文献有差距。
%p1_x0 = [-0.2; 0.2; 0.1]; %粒子轨迹漩涡状波动太小太密。不满意。可能是因为初始位置太高，离漩涡涡流太远。
%p2_x0 = [0.8; 0.15; 0]; %粒子轨迹向上跑。不满意。可能是因为离精子太近（在鞭毛波里面）
%p1_x0 = [-0.2; 0.2; 0]; %粒子轨迹漩涡状波动仍然太小太密，且XY平面上靠近精子的轨迹趋势显弱。不满意。可能是因为离精子太远？
%p2_x0 = [0.8; 0.2; 0];  %粒子轨迹仍然向上跑。不满意。可能是因为在精子Z=0鞭毛波平面内，受向上的漩涡影响。
%
if 1 % Particles near distal tail.
    %x_temp=0.8:0.1:1.0; y_temp=0.18:0.01:0.22; z_temp=-0.1:0.1:0.1;
    x_temp=0.9; y_temp=0.18; z_temp=0;
else % Particles near head
    %x_temp=-0.3:0.1:-0.1; y_temp=0.2:0.1:0.5; z_temp=-0.1:0.1:0.1;
    x_temp=-0.3; y_temp=0.5; z_temp=0.1;
end
nx = length(x_temp); ny = length(y_temp); nz = length(z_temp); 
[X_temp,Y_temp,Z_temp] = ndgrid(x_temp,y_temp,z_temp);  %X_temp: nx*ny*nz 
coords = [X_temp(:); Y_temp(:); Z_temp(:) ];    

    
load XNodes_TetheredSperm_alpha0.4_V2.mat;
dt = tRange(2)-tRange(1);
nt_WF = find(tRange==2*pi)-1;

nt_trace = 100*(size(xs,2)-1);   %100BC
%p1_x = nan(3,nt_trace);
%p2_x = nan(3,nt_trace);
%p1_x(:,1) = p1_x0;
%p2_x(:,1) = p2_x0;
np = nx*ny*nz;  %number of the total particles.
particle_trace = nan(3*np,nt_trace+1);
particle_trace(:,1) = coords;


% numerical parameters
    epsilon=0.25/45;
    domain='i';
    blockSize=0.2;

    for i_nt=1:nt_trace
        i_nt
        
        nt_WF_temp = mod(i_nt,nt_WF);
        if nt_WF_temp==0
            nt_WF_temp = nt_WF;
        end
        xs_temp = xs(:,nt_WF_temp);
        Xs_temp = Xs(:,nt_WF_temp);
        f_temp = f(:,nt_WF_temp);
        
        % Assemble allocation points: sperm and boundary
        x=MergeVectorGrids(xs_temp,xb); 
        X=MergeVectorGrids(Xs_temp,Xb); 
        
        for i_particle = 1:np
            i_particle_coords = [particle_trace(i_particle,i_nt); particle_trace(np+i_particle,i_nt); particle_trace(2*np+i_particle,i_nt)];
                  
            u_p_temp=EvaluateVelocityFromForce(i_particle_coords,X,x,f_temp,epsilon,domain,blockSize);
            
            particle_trace(i_particle,i_nt+1) = particle_trace(i_particle,i_nt) + dt*u_p_temp(1); % x coordinate
            particle_trace(np+i_particle,i_nt+1) = particle_trace(np+i_particle,i_nt) + dt*u_p_temp(2);  % y coordinate
            particle_trace(2*np+i_particle,i_nt+1) = particle_trace(2*np+i_particle,i_nt) + dt*u_p_temp(3);  % z coordinate                   
        end
        
        
        % Evaluate the particle velocity and displacement.
        %u_p1=EvaluateVelocityFromForce(p1_x0,X,x,f_temp,epsilon,domain,blockSize);
        %p1_x(:,i_nt+1) = p1_x(:,i_nt) + dt*u_p1(:);
        %u_p2=EvaluateVelocityFromForce(p2_x0,X,x,f_temp,epsilon,domain,blockSize);
        %p2_x(:,i_nt+1) = p2_x(:,i_nt) + dt*u_p2(:);
    end    
   

save('XNodes_TetheredSperm_alpha0.4_ParticleTrace_Tail_selected.mat','particle_trace','dt');
%save('XNodes_TetheredSperm_alpha0.4_ParticleTrace_Tail.mat','particle_trace','dt');

