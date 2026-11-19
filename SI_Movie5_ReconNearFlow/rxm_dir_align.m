function [xalign,yalign,zalign]=rxm_dir_align(dir_raw,dir_new,x0,y0,z0)
% rotate the raw trajectory ([x0,y0,z0]) in the way that the direction 'dir_rot' is
% aligned with axis X 

% axis-angle representation - coord (vec) (first 3 entries axis)+ angle of ratotion around that axis(for last)
%axisrot = vrrotvec(dir_raw,dir_new); % find rotation axis and angle around which vector a becomes b
%rotmat = vrrotvec2mat(axisrot); % revert axis-angle vector into a rotation matrix
% I need to run this code using 2025 version of MATLAB, so some commands
% from 2017 version's MATLAB need modification, as they are not included
% in the new version, such as 'vrrotvec' and 'vrrotvec2mat'. Below are
% their alternatives.
axisrot = vrrotvec_custom(dir_raw,dir_new);     
rotmat = vrrotvec2mat_custom(axisrot);




% extended rotating matrix         
    on = ones(size(x0));
rotext = {rotmat(1,1)*on,rotmat(1,2)*on,rotmat(1,3)*on;...
          rotmat(2,1)*on,rotmat(2,2)*on,rotmat(2,3)*on;...
          rotmat(3,1)*on,rotmat(3,2)*on,rotmat(3,3)*on};  
      
      
% rotating raw data by the way that makes 'dir_rot' alligned with axis X
X0_R = cell_MatrixProduct(rotext,{x0;y0;z0});
%new trajectory in the aligned coordinate system
        xalign = X0_R{1,1};
        yalign = X0_R{2,1};
        zalign = X0_R{3,1};

