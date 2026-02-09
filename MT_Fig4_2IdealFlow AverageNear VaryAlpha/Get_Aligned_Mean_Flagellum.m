function [Xtail_mean_align, Ytail_mean_align, Ztail_mean_align] = Get_Aligned_Mean_Flagellum(Xtail_0, Ytail_0, Ztail_0, dir_align)

    Ns = size(Xtail_0,1);
    nt_XNodes = size(Xtail_0,2);
    
    % Spatial smooth of the flagellar shape, to remove the noise arising from the raw experimental data being reconstructed.
    ppx = 0.999; % 'ppx=1' means interpolation without cubic spline.
    ppy = 0.999;
    ppz = 0.999;  
    for i_nt = 1:nt_XNodes
        s_temp = arclength(Xtail_0(:,i_nt), Ytail_0(:,i_nt), Ztail_0(:,i_nt));
        Xtail_0(:,i_nt) = fnval(csaps(s_temp,Xtail_0(:,i_nt),ppx),s_temp);
        Ytail_0(:,i_nt) = fnval(csaps(s_temp,Ytail_0(:,i_nt),ppy),s_temp);
        Ztail_0(:,i_nt) = fnval(csaps(s_temp,Ztail_0(:,i_nt),ppz),s_temp);
    end    
    
    
    Xtail_1=Xtail_0-repmat(Xtail_0(1,:),Ns,1);
    Ytail_1=Ytail_0-repmat(Ytail_0(1,:),Ns,1);
    Ztail_1=Ztail_0-repmat(Ztail_0(1,:),Ns,1);   
    
    Xtail_mean= mean(Xtail_1,2); %Ns*1 
    Ytail_mean= mean(Ytail_1,2);
    Ztail_mean= mean(Ztail_1,2);
    
    % Align the flagellum with the specific axis.
    dir_new=dir_align;
    dir_raw=[Xtail_mean(1)-Xtail_mean(end),Ytail_mean(1)-Ytail_mean(end),Ztail_mean(1)-Ztail_mean(end)];
    [Xtail_mean_align,Ytail_mean_align,Ztail_mean_align] = rxm_dir_align(dir_raw,dir_new,Xtail_mean,Ytail_mean,Ztail_mean);
    
    
  