function rxm_SchematicDecayDir_Plot

% This function plots the schematic defining the flow decay directions
% around the sperm body.


clc

load XNodes_sp8_NoBoundStokeslet.mat;
Q=Nhh+Ns;   
    
 
    % Plot.
    fig = 1;
    FAlpha = 1; % Transparency parameter.
    
    figure(fig)
    clf; set(gcf,'color','w');
    hold on; axis equal; axis off;     
    view([140 15])
    

for nt0 = [1 23 130 190 240]% 1:70:1+70*4   % Choose a plot moment.
    
%=====================================================================================================
    % Get sperm shape.
    Xtail_0 = Xs(Nhh+1:Q,nt0); 
    Ytail_0 = Xs(Q+Nhh+1:2*Q,nt0); 
    Ztail_0 = Xs(2*Q+Nhh+1:3*Q,nt0); 
    head_tangent_nt0=head_tangent(nt0,:)';
    head_normal_nt0=head_normal(nt0,:)';
    head_binormal_nt0=cross(head_tangent_nt0,head_normal_nt0);
    B_nt0=[head_tangent_nt0,head_normal_nt0,head_binormal_nt0]; 

    % Get sperm body shape.
    Xtail_1 = repmat(Xtail_0,1,10);
    Ytail_1 = repmat(Ytail_0,1,10);
    Ztail_1 = repmat(Ztail_0,1,10);
    dir_align = [1, 0, 0];
    [Xtail, Ytail, Ztail] = Get_Aligned_Mean_Flagellum(Xtail_1, Ytail_1, Ztail_1, dir_align);   %Xtail_align: Ns*1.  
    [xup1,xup2,xup3,xdown1,xdown2,xdown3,Ind_xup_neck] = Get_CF_Aligned_Head(B_nt0);  % xup1: M1*M2.
    Neck =  [xup1(Ind_xup_neck) xup2(Ind_xup_neck) xup3(Ind_xup_neck)]; 
    Xtail = Xtail-Xtail(2) + Neck(1);  
    Ytail = Ytail-Ytail(2) + Neck(2);
    Ztail = Ztail-Ztail(2) + Neck(3);
    
    
    
    
    %% 
    
   
    % Plot sperm.
    PlotTailCylinder_AimFig(Xtail,Ytail,Ztail,FAlpha,fig)
    s1= surf(xup1,xup2,xup3);  
    s2= surf(xdown1,xdown2,xdown3);
    cc = nan(size(s1.CData,1),size(s1.CData,2),3); cc(:,:,1)=177/255; cc(:,:,2)=24/255; cc(:,:,3)=45/255;s1.CData=cc; 
    cc = nan(size(s1.CData,1),size(s1.CData,2),3); cc(:,:,1)=218/255; cc(:,:,2)=207/255; cc(:,:,3)=168/255;s2.CData=cc; 
    s1.EdgeColor='none' ;
    s2.EdgeColor='none' ;
    s1.FaceAlpha = FAlpha;
    s2.FaceAlpha = FAlpha;
             
   
end



% Plot Cartesian axes.
    O = [-0.45 0 0]; % for PCA modes.
    Cx=[1; 0; 0]; 
    Cy=[0; 1; 0]; 
    Cz=[0; 0; 1]; 
    C=[Cx Cy Cz];
    Plot_Cartesian_Axes(fig,O,C)  
    camlight right
    lighting gouraud    
     
% Plot the six direcitons.
    plot3([-0.45 -0.8],[0 0],[0 0], ':', 'color',[139 0 139]/255,'linewidth',10);  %purple
    plot3([-0.25 0.15],[0 0],[0 0], '-', 'color',[139 0 139]/255,'linewidth',10);  %purple
    plot3([-0.45 -0.45],[0 0],[0 -0.17], ':', 'color',[0 128 0]/255,'linewidth',10);  %green
    plot3([-0.45 -0.45],[0 0],[0.2 0.35], '-', 'color',[0 128 0]/255,'linewidth',10);  %green
    plot3([-0.45 -0.45],[-0.35 0],[0 0], ':', 'color',[255 215 0]/255,'linewidth',10);  %yellow/ golden
    plot3([-0.45 -0.45],[0.2 0.45],[0 0], '-', 'color',[255 215 0]/255,'linewidth',10);  %yellow/ golden
    
    
    
l1=light;
l1.Color = [1 1 1];
l1.Style = 'local';  %'infinite';  %
l1.Position = [0 0.05 1];
l2=light;
l2.Color = [1 1 1];
l2.Style = 'local';  %'infinite';  %
l2.Position = [-0.3 0.1 0.1];
l3=light;
l3.Color = [1 1 1];
l3.Style = 'local';  %'infinite';  %
l3.Position = [-0.6 0.2 0.1];
lighting gouraud   %phong  %flat  %
material metal  %shiny     
