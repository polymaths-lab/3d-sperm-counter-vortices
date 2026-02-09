function [Xtail_mean_align, Ytail_mean_align, Ztail_mean_align] = Get_Aligned_Mean_Flagellum(Xtail_0, Ytail_0, Ztail_0, dir_align)

    Ns = size(Xtail_0,1);
    %nt_XNodes = size(Xtail_0,2);
    
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
    
    
  