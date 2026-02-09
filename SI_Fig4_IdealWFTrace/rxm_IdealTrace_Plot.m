function rxm_IdealTrace_Plot

clc

load XNodes_Lab_alpha0.0.mat;
x_alpha00 = Xs(Nhh+Ns/2,:);
y_alpha00 = Xs(2*Nhh+Ns+Ns/2,:);
z_alpha00 = Xs(3*Nhh+2*Ns+Ns/2,:);
tRange_p = tRange./(2*pi);
% translational velocity data in the lab frame of reference, units:
% flagellum length/ beat cycle.
U_alpha00=Function_derivative_of_time(x_alpha00,tRange_p,'time');
V_alpha00=Function_derivative_of_time(y_alpha00,tRange_p,'time');
W_alpha00=Function_derivative_of_time(z_alpha00,tRange_p,'time');
vel_alpha00 = (U_alpha00.*U_alpha00+V_alpha00.*V_alpha00+W_alpha00.*W_alpha00).^0.5;
y_alpha00(end)=NaN;


load XNodes_Lab_alpha0.2.mat;
x_alpha02 = Xs(Nhh+Ns/2,:);
y_alpha02 = Xs(2*Nhh+Ns+Ns/2,:);
z_alpha02 = Xs(3*Nhh+2*Ns+Ns/2,:);
% translational velocity data in the lab frame of reference, units:
% flagellum length/ beat cycle.
U_alpha02=Function_derivative_of_time(x_alpha02,tRange_p,'time');
V_alpha02=Function_derivative_of_time(y_alpha02,tRange_p,'time');
W_alpha02=Function_derivative_of_time(z_alpha02,tRange_p,'time');
vel_alpha02 = (U_alpha02.*U_alpha02+V_alpha02.*V_alpha02+W_alpha02.*W_alpha02).^0.5;
y_alpha02(end)=NaN;


load XNodes_Lab_alpha0.4.mat;
x_alpha04 = Xs(Nhh+Ns/2,:);
y_alpha04 = Xs(2*Nhh+Ns+Ns/2,:);
z_alpha04 = Xs(3*Nhh+2*Ns+Ns/2,:);
% translational velocity data in the lab frame of reference, units:
% flagellum length/ beat cycle.
U_alpha04=Function_derivative_of_time(x_alpha04,tRange_p,'time');
V_alpha04=Function_derivative_of_time(y_alpha04,tRange_p,'time');
W_alpha04=Function_derivative_of_time(z_alpha04,tRange_p,'time');
vel_alpha04 = (U_alpha04.*U_alpha04+V_alpha04.*V_alpha04+W_alpha04.*W_alpha04).^0.5;
y_alpha04(end)=NaN;


load XNodes_Lab_alpha0.6.mat;
x_alpha06 = Xs(Nhh+Ns/2,:);
y_alpha06 = Xs(2*Nhh+Ns+Ns/2,:);
z_alpha06 = Xs(3*Nhh+2*Ns+Ns/2,:);
% translational velocity data in the lab frame of reference, units:
% flagellum length/ beat cycle.
U_alpha06=Function_derivative_of_time(x_alpha06,tRange_p,'time');
V_alpha06=Function_derivative_of_time(y_alpha06,tRange_p,'time');
W_alpha06=Function_derivative_of_time(z_alpha06,tRange_p,'time');
vel_alpha06 = (U_alpha06.*U_alpha06+V_alpha06.*V_alpha06+W_alpha06.*W_alpha06).^0.5;
y_alpha06(end)=NaN;


load XNodes_Lab_alpha0.8.mat;
x_alpha08 = Xs(Nhh+Ns/2,:);
y_alpha08 = Xs(2*Nhh+Ns+Ns/2,:);
z_alpha08 = Xs(3*Nhh+2*Ns+Ns/2,:);
% translational velocity data in the lab frame of reference, units:
% flagellum length/ beat cycle.
U_alpha08=Function_derivative_of_time(x_alpha08,tRange_p,'time');
V_alpha08=Function_derivative_of_time(y_alpha08,tRange_p,'time');
W_alpha08=Function_derivative_of_time(z_alpha08,tRange_p,'time');
vel_alpha08 = (U_alpha08.*U_alpha08+V_alpha08.*V_alpha08+W_alpha08.*W_alpha08).^0.5;
y_alpha08(end)=NaN;


load XNodes_Lab_alpha1.mat;
x_alpha10 = Xs(Nhh+Ns/2,:);
y_alpha10 = Xs(2*Nhh+Ns+Ns/2,:);
z_alpha10 = Xs(3*Nhh+2*Ns+Ns/2,:);
% translational velocity data in the lab frame of reference, units:
% flagellum length/ beat cycle.
U_alpha10=Function_derivative_of_time(x_alpha10,tRange_p,'time');
V_alpha10=Function_derivative_of_time(y_alpha10,tRange_p,'time');
W_alpha10=Function_derivative_of_time(z_alpha10,tRange_p,'time');
vel_alpha10 = (U_alpha10.*U_alpha10+V_alpha10.*V_alpha10+W_alpha10.*W_alpha10).^0.5;
y_alpha10(end)=NaN;


%% Plot.

fig = 1;

figure(fig)
clf; hold on; axis equal; set(gcf,'color','w')
view([-50 25]); box on; axis off;
XL=[-0.7 0.7]; YL=[-1.8 0.17]; ZL=[0 0.4]; 
xlim([XL(1),XL(2)]);ylim([YL(1),YL(2)]); zlim([ZL(1),ZL(2)]); 

% Lab trace of the specified flagellar point, with velocity values color-encoded.
patch(x_alpha00,y_alpha00+0.05,z_alpha00,vel_alpha00,'EdgeColor','interp','linewidth',3) 
patch(x_alpha02,y_alpha02-0.3,z_alpha02,vel_alpha02,'EdgeColor','interp','linewidth',3) 
patch(x_alpha04,y_alpha04-0.65,z_alpha04,vel_alpha04,'EdgeColor','interp','linewidth',3) 
patch(x_alpha06,y_alpha06-1,z_alpha06,vel_alpha06,'EdgeColor','interp','linewidth',3) 
patch(x_alpha08,y_alpha08-1.32,z_alpha08,vel_alpha08,'EdgeColor','interp','linewidth',3) 
patch(x_alpha10,y_alpha10-1.57,z_alpha10,vel_alpha10,'EdgeColor','interp','linewidth',3) 

% Mark the starting point.
scatter3(x_alpha00(1),y_alpha00(1)+0.05,z_alpha00(1),200,'k','filled')
scatter3(x_alpha02(1),y_alpha02(1)-0.3,z_alpha02(1),200,'k','filled')
scatter3(x_alpha04(1),y_alpha04(1)-0.65,z_alpha04(1),200,'k','filled')
scatter3(x_alpha06(1),y_alpha06(1)-1,z_alpha06(1),200,'k','filled')
scatter3(x_alpha08(1),y_alpha08(1)-1.32,z_alpha08(1),200,'k','filled')
scatter3(x_alpha10(1),y_alpha10(1)-1.57,z_alpha10(1),200,'k','filled')


colormap(parula)
c= colorbar; caxis([0 0.6]); 
set(c,'position',[0.9,0.33,0.007,0.4]);
c.Label.String = 'Trace velocity'; 
c.Ticks=[0 0.6]; c.TickLabels={'0','0.6'}; c.FontName='times'; c.TickLabelInterpreter='latex'; c.FontSize=32;

return
% Plot Cartesian axes.
O=[-0.4 -1.5 0.03];
Cx=[-1; 0; 0]; 
Cy=[0; -1; 0]; 
Cz=[0; 0; 1]; 
C=[Cx Cy Cz];
Plot_Cartesian_Axes(fig,O,C)  
  
camlight right
lighting gouraud
   