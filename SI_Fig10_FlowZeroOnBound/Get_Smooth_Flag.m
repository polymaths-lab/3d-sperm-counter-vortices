function [Xtail_s,Ytail_s,Ztail_s] = Get_Smooth_Flag(Xtail,Ytail,Ztail)


    % Spatial smooth of the flagellar shape, to remove the noise arising from the raw experimental data being reconstructed.
    ppx = 0.999; % 'ppx=1' means interpolation without cubic spline.
    ppy = 0.999;
    ppz = 0.999;  
    
    s_temp = arclength(Xtail,Ytail,Ztail); 
    
    Xtail_s = fnval(csaps(s_temp,Xtail,ppx),s_temp);
    Ytail_s = fnval(csaps(s_temp,Ytail,ppy),s_temp);
    Ztail_s = fnval(csaps(s_temp,Ztail,ppy),s_temp);
    
