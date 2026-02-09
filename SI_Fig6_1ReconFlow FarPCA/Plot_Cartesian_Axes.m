function Plot_Cartesian_Axes(fig,O,C)


    figure(fig)
    hold on;
    
    scatter3(O(1),O(2),O(3),20,'k','filled')
    
    % Vector: line segment.
    eff_x=0.3; eff_y=0.3; eff_z=0.3; 
    quiver3(O(1),O(2),O(3),eff_x*C(1,1),eff_x*C(2,1),eff_x*C(3,1),'k','linewidth',4); % arrow indicating X axis.
    quiver3(O(1),O(2),O(3),eff_y*C(1,2),eff_y*C(2,2),eff_y*C(3,2),'k','linewidth',4); % arrow indicating Y axis.
    quiver3(O(1),O(2),O(3),eff_z*C(1,3),eff_z*C(2,3),eff_z*C(3,3),'k','linewidth',4); % arrow indicating Z axis.
    
    % Vector: arrowhead.
    eff_x=eff_x*0.8; eff_y=eff_y*0.8; eff_z=eff_z*0.8; 
    sz = 0.05; 
    [cone_X_x, cone_X_y, cone_X_z]=meshgrid(-3:3:3, -3:3:3, -1:1:1);  %this may need modification for different scales of figures.
    cone_X_u=repmat(C(1,1),3,3,3);   
    cone_X_v=repmat(C(2,1),3,3,3); 
    cone_X_w=repmat(C(3,1),3,3,3); 
    cone_X_sxyz=[O(1)+eff_x*C(1,1),O(2)+eff_x*C(2,1),O(3)+eff_x*C(3,1)];
    hcone_X = coneplot(cone_X_x, cone_X_y, cone_X_z, cone_X_u, cone_X_v, cone_X_w,cone_X_sxyz(1),cone_X_sxyz(2),cone_X_sxyz(3),sz);
    hcone_X.FaceColor = [0 0 0];
    hcone_X.EdgeColor = 'none';
    hcone_X.DiffuseStrength = 0.8;  
    cone_Y_u=repmat(C(1,2),3,3,3);   
    cone_Y_v=repmat(C(2,2),3,3,3); 
    cone_Y_w=repmat(C(3,2),3,3,3); 
    cone_Y_sxyz=[O(1)+eff_y*C(1,2),O(2)+eff_y*C(2,2),O(3)+eff_y*C(3,2)];
    hcone_Y = coneplot(cone_X_x, cone_X_y, cone_X_z, cone_Y_u, cone_Y_v, cone_Y_w,cone_Y_sxyz(1),cone_Y_sxyz(2),cone_Y_sxyz(3),sz);
    hcone_Y.FaceColor = [0 0 0];
    hcone_Y.EdgeColor = 'none';
    hcone_Y.DiffuseStrength = 0.8;  
    cone_Z_u=repmat(C(1,3),3,3,3);   
    cone_Z_v=repmat(C(2,3),3,3,3); 
    cone_Z_w=repmat(C(3,3),3,3,3); 
    cone_Z_sxyz=[O(1)+eff_z*C(1,3),O(2)+eff_z*C(2,3),O(3)+eff_z*C(3,3)];
    hcone_Z = coneplot(cone_X_x, cone_X_y, cone_X_z, cone_Z_u, cone_Z_v, cone_Z_w,cone_Z_sxyz(1),cone_Z_sxyz(2),cone_Z_sxyz(3),sz);
    hcone_Z.FaceColor = [0 0 0];
    hcone_Z.EdgeColor = 'none';
    hcone_Z.DiffuseStrength = 0.8;
    
    