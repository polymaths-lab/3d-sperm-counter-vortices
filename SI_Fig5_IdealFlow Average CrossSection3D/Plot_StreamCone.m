function Plot_StreamCone(fig,Xg,Yg,Zg,Ug_mean,Vg_mean,Wg_mean,sx_cone,sy_cone,sz_cone,scal)

    figure(fig)
    hold on
    
    hcone = coneplot(Xg,Yg,Zg,Ug_mean,Vg_mean,Wg_mean,sx_cone,sy_cone,sz_cone,scal);
    hcone.FaceColor = 'r';
    hcone.EdgeColor = 'none';
    hcone.DiffuseStrength = 0.8;
    
    
    
    