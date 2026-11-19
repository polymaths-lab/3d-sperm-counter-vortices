function [Xhead_up,Yhead_up,Zhead_up, Xhead_down,Yhead_down,Zhead_down] = Generate_Human_Head_UpDownSides



 % creating mesh for sperm head provide by Hermes
    Ndiv = 60;
    radius_plot = 1.2/56;%normalized
    head_type = 1;    
    [Xhead,Yhead,Zhead] = sperm_head3D_shape(radius_plot, head_type, Ndiv);%size(Zhead)=[61 61]
    scal=1.4; %1.7;%scal=80;%scal=56;  %This 'scal' parameter can be adjusted according to visualization requirement.
    Xhead=Xhead*scal;Yhead=Yhead*scal;Zhead=Zhead*scal;
    M1=size(Xhead,1);M2=size(Xhead,2);M=M1*M2;
    [~,Ind_Neck] = max(Xhead(:)); Xhead_temp = Xhead(:); Xhead_neck = Xhead_temp(Ind_Neck);
   
    % head upside
    [zup_row,zup_col]=find(Zhead>=0); %length(zup_row)=1891
    Zhead_up=nan(size(Zhead));Xhead_up=nan(size(Xhead));Yhead_up=nan(size(Yhead));
    for i = 1:length(zup_row)
        Zhead_up(zup_row(i),zup_col(i))=Zhead(zup_row(i),zup_col(i));
        Xhead_up(zup_row(i),zup_col(i))=Xhead(zup_row(i),zup_col(i));
        Yhead_up(zup_row(i),zup_col(i))=Yhead(zup_row(i),zup_col(i));
    end
    %Xhead_up = Xhead_up-Xhead_neck;
    %xup=[reshape(Xhead_up,M,1);reshape(Yhead_up,M,1);reshape(Zhead_up,M,1)]; 
    
    % head downside
    [zdown_row,zdown_col]=find(Zhead<=0);%length(zdown_row)=1891
    Zhead_down=nan(size(Zhead));Xhead_down=nan(size(Xhead));Yhead_down=nan(size(Yhead));
    for i = 1:length(zdown_row)
        Zhead_down(zdown_row(i),zdown_col(i))=Zhead(zdown_row(i),zdown_col(i));
        Xhead_down(zdown_row(i),zdown_col(i))=Xhead(zdown_row(i),zdown_col(i));
        Yhead_down(zdown_row(i),zdown_col(i))=Yhead(zdown_row(i),zdown_col(i));
    end
    %Xhead_down = Xhead_down-Xhead_neck;
    %xdown=[reshape(Xhead_down,M,1);reshape(Yhead_down,M,1);reshape(Zhead_down,M,1)]; 

    
    
    
 %%   
    
if 0    
% Sperm head:up+down£¬different colors,surf
    model.a1=2.0/45;
    model.a2=1.6/45;
    model.a3=1.0/45;
    nth=20;%10;
    nphi=40;%20;
    a=1;
    %M=nth*nphi;
    
    %head: up side
    [xup1,xup2,xup3]=GenerateSphereSurfaceForVisualisation_up(nth,nphi,a);
    %[M1 M2]=size(xup1);
    xup1=xup1*model.a1;
    xup2=xup2*model.a2;
    xup3=xup3*model.a3;
    
    %head: down side
    [xdown1,xdown2,xdown3]=GenerateSphereSurfaceForVisualisation_down(nth,nphi,a);
    xdown1=xdown1*model.a1;
    xdown2=xdown2*model.a2;
    xdown3=xdown3*model.a3;
end
    