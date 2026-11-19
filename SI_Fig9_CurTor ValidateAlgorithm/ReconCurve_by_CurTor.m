function [s,T,N,B,X] = ReconCurve_by_CurTor(curvature,torsion,T0,N0,X0,sspan)
%% Compute curves by defined curvature and torsion values

% INPUT:
% curvature,torsion: 1*ns
% T0: 1*3. initial tangent vector
% N0: 1*3. initial normal vector
% X0: 1*3. start point of the curve
% sspan: 1*ns or 1*2. integrating range

% OUTPUT:
% s: ns*1. integral interval points, each corresponding to a local Frenet basis [T,N,B]
% T: ns*3. tangent vector
% N: ns*3. normal vector
% B: ns*3. binormal vector
% X: ns*3. spatial coordinates along the target curve 


B0 = cross(T0,N0);
w0 = [T0,N0,B0];
[s,w] = ode45(@(s,w) ODE_SolveFrenetEquations(curvature,torsion,w,s,sspan), sspan, w0);

T = w(:,1:3);
N = w(:,4:6);
B = w(:,7:9);

X_temp = 0.5* cumsum(repmat(diff(s),1,3) .* (T(2:end,:) + T(1:end-1,:)));
ns = size (X_temp,1)+1;
X = nan(ns,3);
X(1,:) = X0;
X(2:end,:) = X_temp+repmat(X0,ns-1,1);


