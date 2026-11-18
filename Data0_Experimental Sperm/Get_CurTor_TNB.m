function [kappa, tau, tt, nn, bb] = Get_CurTor_TNB(x,B0)
%% Measure the curvature, torsion, and local Frenet frames, of the object curve.

% INPUT:
% x: ns*3, spatial coordinates of points along the curve.

ns = size(x,1);
arc = arclength(x(:,1), x(:,2), x(:,3));  % ns*1.    

tt = nan(ns,3);  % tangent vector (unit). valid entries: 1 --> ns-1.
nn = nan(ns,3);  % normal vector (unit). valid entries: 2 --> ns-1.
bb = nan(ns,3);  % binormal vector (unit). valid entries: 2 --> ns-1.
kappa = nan(ns,1);  % curvaure. valid entries: 2 --> ns-1.
tau = nan(ns,1);    % torsion. valid entries: 3 --> ns-1.


% Tangent vector.
drds = diff(x)./diff(arc);  % (ns-1)*3
for i_ns =1:ns-1
    tt(i_ns,:) = drds(i_ns,:)./norm(drds(i_ns,:));
end


% Binormal vector.
for i_ns = 2:ns-1
    bb_temp(i_ns,:) = cross(tt(i_ns-1,:),tt(i_ns,:));
    bb (i_ns,:) = bb_temp(i_ns,:)/norm(bb_temp(i_ns,:));
end
% Sign: specify manually, choose a trihedron of the initial 'bb' pointing along +z direction.
if B0*bb(2,:)'<0  
    bb = -bb;
end
% Enforce frame continuity. Sometimes B flips sign randomly.
for i_ns = 3:ns-1
    if dot(bb(i_ns,:), bb(i_ns-1,:)) < 0
        bb(i_ns,:) = -bb(i_ns,:);
    end
end


% Normal vector.
for i_ns = 2:ns-1
    nn_temp(i_ns,:) = cross(bb(i_ns,:),tt(i_ns,:));
    nn(i_ns,:) = nn_temp(i_ns,:)/norm(nn_temp(i_ns,:));
end


% Curvature.
dr = diff(x); % (ns-1)*3.
for i_ns = 1:ns-2
    dr_1 = dr(i_ns,:);
    dr_2 = dr(i_ns+1,:);
    dr_temp = (norm(dr_1)+norm(dr_2))/2;    
    theta_sin = norm(cross(dr_1,dr_2));
    theta_cos = dot(dr_1,dr_2);
    theta = atan2(theta_sin,theta_cos);  % angles between consective tangent vectors  
    % from Hermes: To convert an integrated quantity back to a pointwise one, we simply divide by the length, |D|, of
    %the domain of integration. In the discrete case, we define Di as the Voronoi region associated to each vertex, having    
    kappa(i_ns+1) = abs(2*tan(theta/2)/dr_temp);
    
    tt_1 = tt(i_ns,:);
    tt_2 = tt(i_ns+1,:);
    dtt = tt_2-tt_1;
    sign_kappa = sign(dot(dtt,nn(i_ns+1,:)));
    kappa(i_ns+1) = kappa(i_ns+1)*sign_kappa;
end


% Torsion. 
for i_ns= 2:ns-2  
    bb_1 = bb(i_ns,:);
    bb_2 = bb(i_ns+1,:);
    phi_sin = norm(cross(bb_1,bb_2));
    phi_cos = dot(bb_1,bb_2);
    phi = atan2(phi_sin,phi_cos);
    dr_1 = dr(i_ns,:);
    dr_2 = dr(i_ns+1,:);
    dr_temp = (norm(dr_1)+norm(dr_2))/2;    
    tau(i_ns+1) = abs(2*tan(phi/2)/dr_temp);
    
    dbb = bb_2-bb_1;
    sign_tau = -1*sign(dot(dbb,nn(i_ns+1,:)));
    tau(i_ns+1) = tau(i_ns+1)*sign_tau;
end



