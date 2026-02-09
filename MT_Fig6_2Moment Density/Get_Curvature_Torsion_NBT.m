function [kappa, tau, tt, nn, bb] = Get_Curvature_Torsion_NBT(x,B0)
%% Measure the curvature and torsion, as well as the instantaneous Frenet frames, of the object curve
%
% INPUT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% t ... nt*1 vector, temporal or
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% spatial integaral
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% points%%%%%%%%%%%%%%%%%%%%%%%%%
% x ... nt*3 array, coordinates of points along the curve over time or space 
% OUTPUT:
% kappa ... 1*(nt-2) vector, the measured curvature, with sign
% tau ... 1*(nt-3) vector, the measured torsion, with sign
% tt ... (nt-1)*3 array, measured tangent vectors
% nn ... (nt-2)*3 array, measured normal vectors
% bb ... (nt-2)*3 array, measured binormal vectors
% Actually, the length of all these OUTPUTs is 'nt', but some data may be
% 'nan' considering the fact that the parameters at the start/end of the
% curve cannot be measured.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note that 'Get_Curvature_Torsion_NBT_V2' is basically the same as
% 'Get_Curvature_Torsion_NBT_V1', except for the simplicity and improvement
% of removing the input 't':
% [kappa, tau, tt, nn, bb] = Get_Curvature_Torsion_NBT_V1(x,t,B0)




%nt = length(t);
%T = repmat(t,1,3);
%x_p = diff(x)./diff(T); 
nt = size(x,1);
x_p = diff(x); 

% Measure tangent vectors===================================================
tt_temp = x_p;  % (nt-1)*3 array, measured tangent vectors
tt = nan(nt,3);  
for i_nt =1:length(tt_temp)
    tt(i_nt,:) = tt_temp(i_nt,:)/norm(tt_temp(i_nt,:));
end


% rxm comment:这里代码可以改进。按照TNB定义来测量，应该先确定N，再由T和N确定B。
% Measure binormal vectors=================================================
bb_temp = nan(length(tt_temp)-1,3);  % (nt-2)*3 array, measured binormal vectors
bb = nan(size(tt));  
for i_nt = 1:length(bb_temp)
    bb_temp(i_nt,:) = cross(tt(i_nt,:),tt(i_nt+1,:));
    bb (i_nt+1,:) = bb_temp(i_nt,:)/norm(bb_temp(i_nt,:));
end
%
%
% define 'sign' (for Frenet vectors 'bb' and 'nn', as well as for 'kappa')===========
Sign = nan(1,nt); % the real size of 'sign' should be [1,nt-2], because sign(1)=sign(end)=nan
%
% specify trihedron orientation manually: choose trihedron such that initial bb points along (+z) direction
if B0*bb(2,:)'<0  %[0 1 0]*bb(2,:)'<0 
    Sign(2) = -1;
else
    Sign(2) = 1;
end
% in case there are sign-flip turns to previous curve segment (especially for twisted space curves)
for i_nt = 3:length(bb)-1
    if bb(i_nt-1,:)*bb(i_nt,:)'<0  
        Sign(i_nt) = -1*Sign(i_nt-1);
    else
        Sign(i_nt) = 1*Sign(i_nt-1);
    end
end
%
for i_nt = 2:length(bb)-1
    bb(i_nt,:) = bb(i_nt,:)*Sign(i_nt); 
end
      
            

% Measure normal vectors===================================================
nn_temp = nan(size(bb_temp));  % (nt-2)*3 array, measured normal vectors
nn = nan(size(tt));  % nt*3 array, including some NAN data
for i_nt = 2:length(nn)-1
    nn_temp(i_nt,:) = cross(bb(i_nt,:),tt(i_nt,:));
    nn(i_nt,:) = nn_temp(i_nt,:)/norm(nn_temp(i_nt,:));
end


% Measure curvature 'kappa'================================================
%kappa_temp = nan(1,nt);
kappa = nan(1,nt); % the real size(kappa)=[1,nt-2] as some 'kappa' data is NAN
%
dsi = nan(1,nt); % processed length of curve segments, and the real size(dsi)=[1,nt-2]
% theta = nan(size(dsi)); % angles between consective tangent vectors, and the real size(theta)=[1,nt-2]
Delta_x = diff(x); % (nt-1)*3 array
for i_nt = 1:length(Delta_x)-1
    dsi_1 = Delta_x(i_nt,:);
    dsi_2 = Delta_x(i_nt+1,:);
    dsi(i_nt+1) = (norm(dsi_1)+norm(dsi_2))/2;
    
    theta_Y = norm(cross(dsi_1,dsi_2));
    theta_X = dot(dsi_1,dsi_2);
    theta = atan2(theta_Y,theta_X);  % angles between consective tangent vectors
    % theta(i_nt+1) = atan2(theta_Y,theta_X);
    
    tt_1 = tt(i_nt,:);
    tt_2 = tt(i_nt+1,:);
    delta_tt = tt_2-tt_1;
    %nn(i_nt+1,:)
    sign_kappa = sign(dot(delta_tt,nn(i_nt+1,:)));

% from Hermes: To convert an integrated quantity back to a pointwise one, we simply divide by the length, |D|, of
%the domain of integration. In the discrete case, we define Di as the Voronoi region associated to each vertex, having    
    kappa(i_nt+1) = abs(2*tan(theta/2)/dsi(i_nt+1))*sign_kappa;
end
%
%kappa_temp = 2*tan(theta./2)./dsi;
% sign
%kappa = kappa_temp.*Sign;



% Measure torsion 'tau'====================================================
tau = nan(1,nt); % real size(tau)=[1,nt-3] due to the lack of some data (NAN)
%
% phi = nan(size(dsi)); % angles between consective binormal vectors without signs, and the real size(phi)=[1,nt-3]
for i_nt= 2:length(bb)-2 % the real length of 'bb' is (nt-2), with the start and end data NAN
    bb_1 = bb(i_nt,:);
    bb_2 = bb(i_nt+1,:);
    phi_Y = norm(cross(bb_1,bb_2));
    phi_X = dot(bb_1,bb_2);
    phi = atan2(phi_Y,phi_X);
    %phi_temp = atan2(phi_Y,phi_X);
% sign of the angle 'phi'  
    %delta_bb = cross(bb_1,bb_2);  
    %sign_phi = sign(dot(delta_bb,tt(i_nt+1,:)));
    %phi(i_nt+1) = sign_phi*phi_temp;
    
    delta_bb = bb_2-bb_1;
    sign_tau = -1*sign(dot(delta_bb,nn(i_nt+1,:)));
    
    tau(i_nt+1) = abs(2*tan(phi./2)./dsi(i_nt+1))*sign_tau;
end


