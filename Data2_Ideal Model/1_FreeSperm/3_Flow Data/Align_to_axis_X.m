function [xalign,yalign,zalign]=Align_to_axis_X(dir_rot,x0,y0,z0)

axisrot = vrrotvec(dir_rot,[-1,0,0]);
rotmat = vrrotvec2mat(axisrot); 


% extended rotational matrix:         
    on = ones(size(x0));
rotext = {rotmat(1,1)*on,rotmat(1,2)*on,rotmat(1,3)*on;...
          rotmat(2,1)*on,rotmat(2,2)*on,rotmat(2,3)*on;...
          rotmat(3,1)*on,rotmat(3,2)*on,rotmat(3,3)*on};  
      
      
% rotate raw data to make 'dir_rot' alligned with axis X.
X0_R = cell_MatrixProduct(rotext,{x0;y0;z0});
        xalign = X0_R{1,1};
        yalign = X0_R{2,1};
        zalign = X0_R{3,1};

