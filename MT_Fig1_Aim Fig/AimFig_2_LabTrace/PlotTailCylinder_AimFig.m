function PlotTailCylinder_AimFig_LabTraj(xt,yt,zt,FAlpha,fig)

% Input: flagellar coordinates. size(X)=ns*1.


% Define the cylinder parameters.   
nz = size(xt,1);
r_temp1=linspace(0,-1.4,0.8*nz);
r_temp2=linspace(-1.4,-3,0.2*nz);
r_temp = [r_temp1 r_temp2];
r = 2.^r_temp;%4.^r_temp;
nr=51;
[X,Y,Z] = cylinder(r*0.007,nr-1); %size(X)=nz*nr.



%% Plot the 3D flagellar shape in the form of a cylinder surface.

figure(fig)
hold on; axis equal; 

% 1. Plot the cylinder.
Cylinder = surf(X,Y,Z);


% 2. Rotate the cyclinder.
% Calculate the rotation information of the cylinder.
ax_cylinder = [0 0 1]; % Original cylinder direction.
ax_neck = [xt(2)-xt(1), yt(2)-yt(1), zt(2)-zt(1),];% Flagellar neck vector.
angle_rot = acos(dot(ax_cylinder,ax_neck ) /(norm(ax_cylinder)*norm(ax_neck ))) *180/pi;  %Cylinder rotation angle.
ax_rot = cross(ax_cylinder, ax_neck);  %Cylinder rotation axis.
O_rot = [0 0 0]; %Origin of the cylinder rotation.
% Rotate the cylinder to align with the flagellar neck vector.
rotate(Cylinder, ax_rot, angle_rot, O_rot);


% 3. Translate the cylinder.
xx = Cylinder.XData; x0 = mean(xx,2);  xx = xx - repmat(x0,1,nr) + repmat(xt,1,nr);
yy = Cylinder.YData; y0 = mean(yy,2);  yy = yy - repmat(y0,1,nr) + repmat(yt,1,nr);
zz = Cylinder.ZData; z0 = mean(zz,2);  zz = zz - repmat(z0,1,nr) + repmat(zt,1,nr);
Cylinder.XData = xx;
Cylinder.YData = yy;
Cylinder.ZData = zz;


             
% 4. Prescribe the cylinder outlook, such as color and transparency.
cc = nan(size(Cylinder.CData,1),size(Cylinder.CData,2),3); 
cc(:,:,1)=36/255; cc(:,:,2)=100/255; cc(:,:,3)=171/255; Cylinder.CData=cc;  
Cylinder.FaceAlpha = FAlpha;

shading interp



%% Plot the 2D projected flagellar shape in the form of a cylinder surface.

% 1. Plot the cylinder.
Cylinder = surf(X,Y,Z);


% 2. Rotate the cyclinder.
% Calculate the rotation information of the cylinder.
ax_cylinder = [0 0 1]; % Original cylinder direction.
ax_neck = [xt(2)-xt(1), yt(2)-yt(1), zt(2)-zt(1),];% Flagellar neck vector.
angle_rot = acos(dot(ax_cylinder,ax_neck ) /(norm(ax_cylinder)*norm(ax_neck ))) *180/pi;  %Cylinder rotation angle.
ax_rot = cross(ax_cylinder, ax_neck);  %Cylinder rotation axis.
O_rot = [0 0 0]; %Origin of the cylinder rotation.
% Rotate the cylinder to align with the flagellar neck vector.
rotate(Cylinder, ax_rot, angle_rot, O_rot);


% 3. Translate the cylinder.
xx = Cylinder.XData; x0 = mean(xx,2);  xx = xx - repmat(x0,1,nr) + repmat(xt,1,nr);
yy = Cylinder.YData; y0 = mean(yy,2);  yy = yy - repmat(y0,1,nr) + repmat(yt,1,nr);
zz = Cylinder.ZData; z0 = mean(zz,2);  zz = zz - repmat(z0,1,nr) + repmat(zt,1,nr);
Cylinder.XData = xx;
Cylinder.YData = yy;
Cylinder.ZData = zz;


             
% 4. Prescribe the cylinder outlook, such as color and transparency.
cc = nan(size(Cylinder.CData,1),size(Cylinder.CData,2),3); cc(:,:,1)=0; cc(:,:,2)=0; cc(:,:,3)=0;Cylinder.CData=cc;
Cylinder.FaceAlpha = FAlpha;


% 5. Set Z=0 for projected flagellum.
h = Cylinder;
h.ZData = h.ZData.*0-0.3;


shading interp


