function rxm_Empty_ExpCellData
% This function processes the experimental data to make all cells' data empty, except for the selected one (sp#8).
clc

%%
if 0
load 03_Head_Fixed_Frame_Method_Rotations_Traces_raw_5_Microns_Difraction_head_spin_rxmHeadNBT.mat;
for sp=[1:7 9:30]
    X{sp}=[];  
    Y{sp}=[];  
    Z{sp}=[];  
    head_spin_angle{sp}=[]; 
    sperm_timePointsAnalyzed{sp}=[]; 
    head_binormal{sp}=[]; 
    head_tangent{sp}=[]; 
    head_normal{sp}=[]; 
end
end

%%
if 0
    load FreeSperm_Frequency.mat;
    for sp=[1:7 9:30]
        CF_Freq(sp)=NaN;
        HF_Freq(sp)=NaN;
        HeadSpin_Freq(sp)=NaN;
        NF_Freq(sp)=NaN;
    end
end


%%
clearvars sp
save('new.mat');