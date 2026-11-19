function [Traj_X, Traj_X_new, CA_traj] = Align_traj_with_AP_CA(x0,y0,z0,nspiral)

% This function is to align raw trajectory to the postion where the central
% axis conincides with negative axis X. During this process, we need to
% calculate average path (AP) and central axis (CA) in sequence, and then
% rotate and translate raw paths.


% averaged path, using least-square method.
nt = size(x0,1);
n_points=round(nt/4);
pp1=0.9999;
[x1_temp,y1_temp,z1_temp,~] = fn_interpolate_arclength(x0,y0,z0,n_points,pp1); % to calculate AP, where the smoothing parameter 'pp' varies for planar and 3d raw trajectories.
pp2=1;
[x1,y1,z1,~] = fn_interpolate_arclength(x1_temp,y1_temp,z1_temp,nt,pp2); % to make the AP points denser, where 'pp=1', performing no smoothing function.


% centre axis (of helix), using sliding windows
window= round(nt/nspiral); 
i_x2=0;
for i=[window length(x1)]  
    i_x2= i_x2+1;
    x2_temp(i_x2)=sum(x1((i-window+1):i))/window;  
    y2_temp(i_x2)=sum(y1((i-window+1):i))/window;  
    z2_temp(i_x2)=sum(z1((i-window+1):i))/window;  
end

x2 = linspace(x2_temp(1),x2_temp(end),nt)';
y2 = linspace(y2_temp(1),y2_temp(end),nt)';
z2 = linspace(z2_temp(1),z2_temp(end),nt)';


% align and translate the raw trajectory, as well as the AP, so that the aligned center axis coincides with axis X.
dir_rot=[x2_temp(end)-x2_temp(1),y2_temp(end)-y2_temp(1),z2_temp(end)-z2_temp(1)];
[x0_new_temp,y0_new_temp,z0_new_temp]=Align_to_axis_X(dir_rot,x0,y0,z0); % align with negative axis X.
[x1_new_temp,y1_new_temp,z1_new_temp]=Align_to_axis_X(dir_rot,x1,y1,z1); % align with negative axis X.
[~,y2_new_temp,z2_new_temp]=Align_to_axis_X(dir_rot,x2,y2,z2);
%
x0_new = x0_new_temp-x0_new_temp(1);
y0_new = y0_new_temp-y2_new_temp(1);
z0_new = z0_new_temp-z2_new_temp(1);
%
x1_new = x1_new_temp-x0_new_temp(1);
y1_new = y1_new_temp-y2_new_temp(1);
z1_new = z1_new_temp-z2_new_temp(1);
%
x2_new = linspace(x0_new(1),x0_new(end),nt)';
y2_new = x2_new*0;
z2_new = x2_new*0; 


Traj_X = [x0,y0,z0,x1,y1,z1,x2,y2,z2];
Traj_X_new = [x0_new,y0_new,z0_new,x1_new,y1_new,z1_new,x2_new,y2_new,z2_new];
CA_traj = dir_rot/norm(dir_rot);
