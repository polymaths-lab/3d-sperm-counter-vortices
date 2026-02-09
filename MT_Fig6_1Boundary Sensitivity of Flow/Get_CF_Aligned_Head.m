function [xup1_nt0_align,xup2_nt0_align,xup3_nt0_align, xdown1_nt0_align,xdown2_nt0_align,xdown3_nt0_align, Ind_xup_neck, Ind_xup_front] = Get_CF_Aligned_Head(B_nt0)



% Generate head surface.
[xup1,xup2,xup3, xdown1,xdown2,xdown3] = Generate_Human_Head_UpDownSides;  %xup1: M1*M2.
[~,Ind_xup_neck] = max(xup1(:));
[~,Ind_xup_front] = min(xup1(:)); 
HC = ([xup1(Ind_xup_front) xup2(Ind_xup_front) xup3(Ind_xup_front)] + [xup1(Ind_xup_neck) xup2(Ind_xup_neck) xup3(Ind_xup_neck)])/2;
xup1 = xup1-HC(1);
xdown1 = xdown1-HC(1);


[M1, M2] = size(xup1);
M = M1*M2;
xup=[reshape(xup1,M,1);reshape(xup2,M,1);reshape(xup3,M,1)]; 
xdown=[reshape(xdown1,M,1);reshape(xdown2,M,1);reshape(xdown3,M,1)]; 
  
      

% Align the head.
% head surface rotation and translation.
% upside head
xup_nt0_align=ApplyRotationMatrix(B_nt0,xup); % This is for the lab frame.    
%xup_nt0_align=TranslatePoints(xup_nt0_align,HC);
[xup1_nt0_align,xup2_nt0_align,xup3_nt0_align]=ExtractComponents(xup_nt0_align);
xup1_nt0_align=reshape(xup1_nt0_align,M1,M2);
xup2_nt0_align=reshape(xup2_nt0_align,M1,M2);
xup3_nt0_align=reshape(xup3_nt0_align,M1,M2);
% downside head
xdown_nt0_align=ApplyRotationMatrix(B_nt0,xdown);  % This is for the lab frame.
%xdown_nt0_align=TranslatePoints(xdown_nt0_align,HC);
[xdown1_nt0_align,xdown2_nt0_align,xdown3_nt0_align]=ExtractComponents(xdown_nt0_align);
xdown1_nt0_align=reshape(xdown1_nt0_align,M1,M2);
xdown2_nt0_align=reshape(xdown2_nt0_align,M1,M2);
xdown3_nt0_align=reshape(xdown3_nt0_align,M1,M2);
    
    
    



