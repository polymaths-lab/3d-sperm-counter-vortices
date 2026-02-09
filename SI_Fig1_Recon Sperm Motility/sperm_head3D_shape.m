function [x,y,z] = sperm_head3D_shape(r, head_type, Ndiv)
 
 
 
% Making a unitary sphere
  
[X,Y,Z] = sphere(Ndiv); 

    
  % Spherical head
     if (head_type == 0)  
           x_axis =  r;
           y_axis =  r;
 
           %x = x_axis.*x - x_axis; to plot relative to the centroid
           x = x_axis.*X;
           y = y_axis.*Y;
           z = y_axis.*Z;
           
         
    % a Humanoid sperm head
      elseif (head_type == 1)   
            x_axis =  0.5*4.5/56.0;
            y_axis =  0.5*2.8/56.0;
            z_axis =  y_axis*0.4;
 
 
            Xhead = x_axis*X;
            Yhead = y_axis*Y;
            Zhead = z_axis*Z;
 
            xd = Xhead; 
            yd = Yhead;
            zd = Zhead;
 
 
            zd = zd.*( 0.5.*xd./x_axis + 1.5);
            yd = yd.*(-0.2.*xd./x_axis + 0.8);
            
            
            clear Xhead Yhead Zhead
            x = xd-x_axis/1.1;
            y = yd;
            z = zd;
            
         % Tryp Head head - enlongated ellipsoid
      elseif (head_type == 2)   
   
            x_axis =  2*r;
            y_axis =  r;
            z_axis =  y_axis;
 
 
            Xhead = x_axis*X;
            Yhead = y_axis*Y;
            Zhead = z_axis*Z;
 
            xd = Xhead-x_axis; 
            yd = Yhead;
            zd = Zhead;
 
 
            %zd = zd.*( 0.5.*xd./x_axis + 1.5);
            %yd = yd.*(-0.2.*xd./x_axis + 0.8);
            
            
            clear Xhead Yhead Zhead
            x = xd;
            y = yd;
            z = zd;
          
          
  
            
     end
     
%     head_surf_center = [-0.0365 0 0];     
%     x = x-head_surf_center(1);
%     y = y-head_surf_center(2);
%     z = z-head_surf_center(3);
end