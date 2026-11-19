function rxm_Speed_Data
% This function measures sperm speed data: 
% 3 linear and 3 angular speeds (both exp. and num. sperm);
% cross-correlation between the exp. and num. speeds;
% speed frequencies.

clc

sp = 8;

%% Measure experimental speeds.

load 03_Head_Fixed_Frame_Method_Rotations_Traces_raw_5_Microns_Difraction_head_spin_rxmHeadNBT.mat;
b1_exp = -head_normal{sp};  % normalized. 3*nt.
b2_exp = head_binormal{sp}; 
b3_exp = -head_tangent{sp}; 
nt_exp=size(b3_exp,2); 
dt=1/90; 
t_exp=0:dt:(nt_exp-1)*dt;


%==============================================================================================================
% Measure the 3 angular velocities of the sperm head.
db1_exp_dt=Function_derivative_of_time(b1_exp,t_exp);
db2_exp_dt=Function_derivative_of_time(b2_exp,t_exp);
db3_exp_dt=Function_derivative_of_time(b3_exp,t_exp); %3*nt.
wbody_exp=nan(3,nt_exp);
for i_t=1:nt_exp
    wbody_exp(1,i_t)=dot(  db2_exp_dt(:,i_t),b3_exp(:,i_t)  ); % angular speed around head normal vector, b1.
    wbody_exp(2,i_t)=dot(  db3_exp_dt(:,i_t),b1_exp(:,i_t)  ); % angular speed around head binormal vector, b2.
    wbody_exp(3,i_t)=dot(  db1_exp_dt(:,i_t),b2_exp(:,i_t)  ); % angular speed around head tangent vector, b3.
end


%==============================================================================================================
% Measure the 3 linear velocities. 
% Because we don't have experimental tracking of the head (center) positions, we have to approximate the head rigid body's linear velocities according to the neck point and the head orientation.
load 01_lab_frame_raw_traces_2017_2018.mat;
head_a1 =  0.5*4.5; % semi-axis length of the longitudinal head axis.
X0 = [X{sp}(1,:); Y{sp}(1,:); Z{sp}(1,:)] + b3_exp*head_a1;  % approximate head center coordinates.
Ux=Function_derivative_of_time(X0(1,:),t_exp);
Uy=Function_derivative_of_time(X0(2,:),t_exp);
Uz=Function_derivative_of_time(X0(3,:),t_exp);
ubody_exp=nan(3,nt_exp);
for i_nt=1:nt_exp
    R_temp=[b1_exp(:,i_nt)';b2_exp(:,i_nt)';b3_exp(:,i_nt)'];
    u_temp=R_temp*[Ux(i_nt);Uy(i_nt);Uz(i_nt)];
    ubody_exp(:,i_nt)=u_temp;    
end


clear X Y Z


%% Measure numerical speeds (no boundary).

load XNodes_sp8_NoBoundStokeslet.mat;  
X = Xs(Nhh+1:Nhh+Ns,:); 
Y = Xs(2*Nhh+Ns+1:2*Nhh+2*Ns,:); 
Z = Xs(3*Nhh+2*Ns+1:3*Nhh+3*Ns,:);
ns=size(X,1);
nt_num_NB=size(X,2); 
head_tangent=-head_tangent';
head_normal=-head_normal';
head_binormal=nan(size(head_tangent)); %3*nt.
for i_nt=1:nt_num_NB
    head_binormal(:,i_nt)=cross(head_tangent(:,i_nt),head_normal(:,i_nt));          
end

% Time interpolation to get the reconstructed flagellum positions at the experimentally observed time points.
s_inp = linspace(0,1,ns);
t_dim = t_nd/(2*pi)/Freq_WF;
[T_exp,S_exp] = meshgrid(t_exp,s_inp) ;
[T_num,S_num] = meshgrid(t_dim,s_inp) ; 
X = interp2(T_num,S_num,X,T_exp,S_exp,'spline'); 
Y = interp2(T_num,S_num,Y,T_exp,S_exp,'spline'); 
Z = interp2(T_num,S_num,Z,T_exp,S_exp,'spline');
% Time interpolation to get the reconstructed head orientations at the experimentally observed time points.
[T_exp,S_exp] = meshgrid(t_exp,[1 2 3]) ;
[T_num,S_num] = meshgrid(t_dim,[1 2 3]) ; 
b3_num_NB = interp2(T_num,S_num,head_tangent,T_exp,S_exp,'spline');             
b1_num_NB  = interp2(T_num,S_num,head_normal,T_exp,S_exp,'spline');             
b2_num_NB = interp2(T_num,S_num,head_binormal,T_exp,S_exp,'spline');  



%==============================================================================================================
% Measure the 3 angular velocities of the sperm head (neck flagellar point).
% Assuming the first flagellar point, the neck point, is fixed to the head, so these 3 velocities can also be taken as of the neck point, consistent with the 3 linear velocities of the same neck point.
db1_num_NB_dt = Function_derivative_of_time(b1_num_NB,t_exp); %3*nt.
db2_num_NB_dt = Function_derivative_of_time(b2_num_NB,t_exp);
db3_num_NB_dt = Function_derivative_of_time(b3_num_NB,t_exp);
wbody_num_NB=nan(3,nt_exp);
for i_t=1:nt_exp
    wbody_num_NB(1,i_t)=dot(  db2_num_NB_dt(:,i_t),b3_num_NB(:,i_t)  ); % rotational speed around head normal vector, b1.
    wbody_num_NB(2,i_t)=dot(  db3_num_NB_dt(:,i_t),b1_num_NB(:,i_t)  ); % rotational speed around head binormal vector, b2.
    wbody_num_NB(3,i_t)=dot(  db1_num_NB_dt(:,i_t),b2_num_NB(:,i_t)  ); % rotational speed around head tangent vector, b3.
end



%==============================================================================================================
% Measure the 3 linear velocities of the neck flagellar point.
head_a1=2.0/45; % semi-axis length of the longitudinal head axis. Note that the exp. sperm and virtual sperm have different head geometories.
X0 = [X(1,:); Y(1,:); Z(1,:)] + b3_num_NB*head_a1; % approximate head center coordinates.
Ux=Function_derivative_of_time(X0(1,:),t_exp);
Uy=Function_derivative_of_time(X0(2,:),t_exp);
Uz=Function_derivative_of_time(X0(3,:),t_exp);

load 01_lab_frame_raw_traces_2017_2018_rxmSSppxyz.mat; 
arc_mat = arclength(X{sp},Y{sp},Z{sp});
arc_mean = mean(arc_mat(end,:));

ubody_num_NB = nan(3,nt_exp);
for i_nt=1:nt_exp
    R_temp=[b1_num_NB(:,i_nt)';b2_num_NB(:,i_nt)';b3_num_NB(:,i_nt)'];
    u_temp=R_temp*[Ux(i_nt);Uy(i_nt);Uz(i_nt)];
    ubody_num_NB(:,i_nt)=u_temp*arc_mean;
end

clear X Y Z



%% Measure numerical speeds (with boundary).

load XNodes_sp8_BoundBlakelet_H0.2.mat; 
X = Xs(Nhh+1:Nhh+Ns,:); 
Y = Xs(2*Nhh+Ns+1:2*Nhh+2*Ns,:); 
Z = Xs(3*Nhh+2*Ns+1:3*Nhh+3*Ns,:);
ns=size(X,1);
nt_num_WB=size(X,2); 
head_tangent=-head_tangent';
head_normal=-head_normal';
head_binormal=nan(size(head_tangent)); %3*nt.
for i_nt=1:nt_num_WB
    head_binormal(:,i_nt)=cross(head_tangent(:,i_nt),head_normal(:,i_nt));          
end

% Time interpolation to get the reconstructed flagellum positions at the experimentally observed time points.
s_inp = linspace(0,1,ns);
t_dim = t_nd/(2*pi)/Freq_WF;
[T_exp,S_exp] = meshgrid(t_exp,s_inp) ;
[T_num,S_num] = meshgrid(t_dim,s_inp) ; 
X = interp2(T_num,S_num,X,T_exp,S_exp,'spline'); 
Y = interp2(T_num,S_num,Y,T_exp,S_exp,'spline'); 
Z = interp2(T_num,S_num,Z,T_exp,S_exp,'spline');
% Time interpolation to get the reconstructed head orientations at the experimentally observed time points.
[T_exp,S_exp] = meshgrid(t_exp,[1 2 3]) ;
[T_num,S_num] = meshgrid(t_dim,[1 2 3]) ; 
b3_num_WB = interp2(T_num,S_num,head_tangent,T_exp,S_exp,'spline');             
b1_num_WB  = interp2(T_num,S_num,head_normal,T_exp,S_exp,'spline');             
b2_num_WB = interp2(T_num,S_num,head_binormal,T_exp,S_exp,'spline');  



%==============================================================================================================
% Measure the 3 angular velocities of the sperm head (neck flagellar point).
% Assuming the first flagellar point, the neck point, is fixed to the head, so these 3 velocities can also be taken as of the neck point, consistent with the 3 linear velocities of the same neck point.
db1_num_WB_dt = Function_derivative_of_time(b1_num_WB,t_exp); %3*nt.
db2_num_WB_dt = Function_derivative_of_time(b2_num_WB,t_exp);
db3_num_WB_dt = Function_derivative_of_time(b3_num_WB,t_exp);
wbody_num_WB=nan(3,nt_exp);
for i_t=1:nt_exp
    wbody_num_WB(1,i_t)=dot(  db2_num_WB_dt(:,i_t),b3_num_WB(:,i_t)  ); % rotational speed around head normal vector, b1.
    wbody_num_WB(2,i_t)=dot(  db3_num_WB_dt(:,i_t),b1_num_WB(:,i_t)  ); % rotational speed around head binormal vector, b2.
    wbody_num_WB(3,i_t)=dot(  db1_num_WB_dt(:,i_t),b2_num_WB(:,i_t)  ); % rotational speed around head tangent vector, b3.
end



%==============================================================================================================
% Measure the 3 linear velocities of the neck flagellar point.
head_a1=2.0/45; % semi-axis length of the longitudinal head axis. Note that the exp. sperm and virtual sperm have different head geometories.
X0 = [X(1,:); Y(1,:); Z(1,:)] + b3_num_WB*head_a1; % approximate head center coordinates.
Ux=Function_derivative_of_time(X0(1,:),t_exp);
Uy=Function_derivative_of_time(X0(2,:),t_exp);
Uz=Function_derivative_of_time(X0(3,:),t_exp);

load 01_lab_frame_raw_traces_2017_2018_rxmSSppxyz.mat; 
arc_mat = arclength(X{sp},Y{sp},Z{sp});
arc_mean = mean(arc_mat(end,:));

ubody_num_WB = nan(3,nt_exp);
for i_nt=1:nt_exp
    R_temp=[b1_num_WB(:,i_nt)';b2_num_WB(:,i_nt)';b3_num_WB(:,i_nt)'];
    u_temp=R_temp*[Ux(i_nt);Uy(i_nt);Uz(i_nt)];
    ubody_num_WB(:,i_nt)=u_temp*arc_mean;
end



%% Measure speed frequency.

FS = 90;  % Sampling frequency. dt=1/90.
FD = FS*(0:round(nt_exp/2))/nt_exp; % Frequency domain.
[Fre_w1, ~, ~, PS_w1] = Get_Vector_Calibrated_FFT(wbody_exp(1,:), 2, FS);
[Fre_w2, ~, ~, PS_w2] = Get_Vector_Calibrated_FFT(wbody_exp(2,:), 2, FS);
[Fre_w3, ~, ~, PS_w3] = Get_Vector_Calibrated_FFT(wbody_exp(3,:), 2, FS);
[Fre_u1, ~, ~, PS_u1] = Get_Vector_Calibrated_FFT(ubody_exp(1,:), 2, FS);
[Fre_u2, ~, ~, PS_u2] = Get_Vector_Calibrated_FFT(ubody_exp(2,:), 2, FS);
[Fre_u3, ~, ~, PS_u3] = Get_Vector_Calibrated_FFT(ubody_exp(3,:), 2, FS);
    
PS_wu_exp = [PS_w1; PS_w2; PS_w3; PS_u1; PS_u2; PS_u3];
Fre_wu_exp = [Fre_w1; Fre_w2; Fre_w3; Fre_u1; Fre_u2; Fre_u3];


%% Calculate the cross correlation.

wu_exp = [wbody_exp; ubody_exp]; %6*nt
wu_num_NB = [wbody_num_NB; ubody_num_NB];
wu_num_WB = [wbody_num_WB; ubody_num_WB];

cc_NB = nan(6,2*nt_exp-1);
lags_NB = nan(6,2*nt_exp-1);
cc_WB = nan(6,2*nt_exp-1);
lags_WB = nan(6,2*nt_exp-1);
for i_vel = 1:6
    [cc_temp, lags_temp] = xcorr(wu_exp(i_vel,:)',wu_num_NB(i_vel,:)','coeff');
    cc_NB(i_vel,:) = cc_temp(:)';
    lags_NB(i_vel,:) = lags_temp(:)';
    clear cc_temp lags_temp
    
    [cc_temp, lags_temp] = xcorr(wu_exp(i_vel,:)',wu_num_WB(i_vel,:)','coeff');
    cc_WB(i_vel,:) = cc_temp(:)';
    lags_WB(i_vel,:) = lags_temp(:)';
end



%% Save data

wu_num_B02 = wu_num_NB;
cc_B02 = cc_WB;
lags_B02 = lags_WB;

save('Speed.mat','t_exp','wu_exp','wu_num_NB','wu_num_B02','FS','FD','PS_wu_exp','Fre_wu_exp','cc_NB','lags_NB','cc_B02','lags_B02');