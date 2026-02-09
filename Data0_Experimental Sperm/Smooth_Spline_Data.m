function  Smooth_Spline_Data

% This function is used to smooth the raw data of flagellar beating in
% various frames of reference.


clc

%load 03_Head_Fixed_Frame_Method_Rotations_Traces_raw_5_Microns_Difraction_head_spin.mat;
load 01_lab_frame_raw_traces_2017_2018.mat;
xi = cell(1,30);
yi = cell(1,30);
zi = cell(1,30);
for sp = 8 %1:30

        x_raw{sp} = X{sp};
        y_raw{sp} = Y{sp};
        z_raw{sp} = Z{sp};        
    np_raw(sp)   = size(x_raw{sp},1);
    nt_raw(sp)   = size(x_raw{sp},2);


      ppx = 0.1; % 1= interpolation without cubic spline
      ppy = 0.1;
      ppz = 0.01;
      n_points = 100; 

xi{sp} = zeros(n_points, nt_raw(sp)); 
yi{sp} = xi{sp};
zi{sp} = xi{sp};
ai{sp} = xi{sp};

for j = 1:nt_raw(sp)
      clear xtemp ytemp ztemp s_temp
      xtemp = x_raw{sp}(:,j);
      ytemp = y_raw{sp}(:,j);
      ztemp = z_raw{sp}(:,j);
      
      sold = arclength(xtemp, ytemp, ztemp);
      snew = linspace(0,sold(end),n_points)';

xi{sp}(:,j) = fnval(csaps(sold,xtemp,ppx),snew);
yi{sp}(:,j) = fnval(csaps(sold,ytemp,ppy),snew);
zi{sp}(:,j) = fnval(csaps(sold,ztemp,ppz),snew);

end

% EVALUATE NEW ARCLENGTH
%ai{sp} = arclength(xi{sp},yi{sp},zi{sp});


% For neck and head frames.
if 0
xi{sp} = xi{sp} - repmat(xi{sp}(1,:),n_points,1);   
yi{sp} = yi{sp} - repmat(yi{sp}(1,:),n_points,1);   
zi{sp} = zi{sp} - repmat(zi{sp}(1,:),n_points,1);  
end 

end
clear X Y Z
X = xi;
Y = yi;
Z = zi;
%arcl = ai;



%save('03_Head_Fixed_Frame_Method_Rotations_Traces_raw_5_Microns_Difraction_head_spin_rxmSSppxyz.mat','X','Y','Z');
save('01_lab_frame_raw_traces_2017_2018_rxmSSppxyz.mat','X','Y','Z');

