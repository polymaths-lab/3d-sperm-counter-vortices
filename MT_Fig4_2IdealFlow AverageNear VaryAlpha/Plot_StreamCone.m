function Plot_StreamCone(fig,Xg,Yg,Zg,Ug_mean,Vg_mean,Wg_mean,sx_cone,sy_cone,sz_cone,scal)

    figure(fig)
    hold on
    
    
    sx_cone(2,:,2)=NaN;
    sy_cone(2,:,2)=NaN;
    sz_cone(2,:,2)=NaN;
    
    sx_cone([1 3],:,[1 3])=NaN;
    sy_cone([1 3],:,[1 3])=NaN;
    sz_cone([1 3],:,[1 3])=NaN;
    
    hcone = coneplot(Xg,Yg,Zg,Ug_mean,Vg_mean,Wg_mean,sx_cone,sy_cone,sz_cone,scal);
    hcone.FaceColor = 'r';
    hcone.EdgeColor = 'none';
    hcone.DiffuseStrength = 0.8;
    
    
    
    