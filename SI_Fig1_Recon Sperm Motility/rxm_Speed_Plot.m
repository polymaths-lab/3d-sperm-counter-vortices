function rxm_Speed_Plot

% This function plots: 
% 1. The 6 speeds of the experimental baseline and the two numerical recostructions, respectively. 
% 2. Power spectra of the 6 experimental speeds are provided to determine their frequencies. 
% 3. Cross correlations between the baseline and reconstructed speeds are ploted, too.
% 4. Schematic defining the 6 speeds.


clc

sp = 8;

load Speed_SI.mat;

wu_num_1 = wu_num_BH05;
wu_num_2 = wu_num_BH1;
cc_1 = cc_BH05;
lags_1 = lags_BH05;
cc_2 = cc_BH1;
lags_2 = lags_BH1;


cc_1_max = max(cc_1,[],2);
cc_2_max = max(cc_2,[],2);
cc_1_max = roundn(cc_1_max,-3);
cc_2_max = roundn(cc_2_max,-3);


%% Plot the 6 speeds.

if 1
    
fs=34;
 
figure(1)
clf; 
set(gcf,'color','w','units','normalized','outerposition',[0 0 1 1]) 

% Angular speeds
subplot(3,2,1)
hold on; box on;
ax = gca; ax.TickLabelInterpreter = 'latex'; ax.FontSize = fs;
p1 = plot(t_exp,wu_exp(1,:),'k','linewidth',8);
p2 = plot(t_exp,wu_num_1(1,:),'r','linewidth',4);
p3 = plot(t_exp,wu_num_2(1,:),'g*','linewidth',10);
text(0.2,100,num2str(cc_1_max(1)),'color',[178 34 34]./255,'FontSize',fs,'FontName','time','interpreter','latex') 
text(0.7,100,num2str(cc_2_max(1)),'color',[61 145 64]./255,'FontSize',fs,'FontName','time','interpreter','latex') 
if 1
    lgd = legend([p1 p2 p3],{'$Exp.$','$Rec.(0.5L)$','$Rec.(L)$'},'interpreter','latex'); 
    lgd.ItemTokenSize = [160 40];
    set(lgd,'FontSize',fs,'FontName','time','interpreter','latex');
    legend('boxoff')
end
xlim([0 1]);  ylim([-65 150])
xticks(0:0.2:1); yticks([0 100]);
xticklabels([]); ylabel('$\Omega_1 \; (rad \cdot s^{-1})$','FontSize',fs,'FontName','time','interpreter','latex')


subplot(3,2,3)
hold on; box on;
ax = gca; ax.TickLabelInterpreter = 'latex'; ax.FontSize = fs;
plot(t_exp,wu_exp(2,:),'k','linewidth',4);
plot(t_exp,wu_num_1(2,:),'r','linewidth',2);
plot(t_exp,wu_num_2(2,:),'g*','linewidth',1);
text(0.2,100,num2str(cc_1_max(2)),'color',[178 34 34]./255,'FontSize',fs,'FontName','time','interpreter','latex') 
text(0.7,100,num2str(cc_2_max(2)),'color',[61 145 64]./255,'FontSize',fs,'FontName','time','interpreter','latex') 
xlim([0 1]);  ylim([-50 150])
xticks(0:0.2:1); yticks([0 100]);
xticklabels([]); ylabel('$\Omega_2 \; (rad \cdot s^{-1})$','FontSize',fs,'FontName','time','interpreter','latex')


subplot(3,2,5)
hold on; box on;
ax = gca; ax.TickLabelInterpreter = 'latex'; ax.FontSize = fs;
plot(t_exp,wu_exp(3,:),'k','linewidth',4);
plot(t_exp,wu_num_1(3,:),'r','linewidth',2);
plot(t_exp,wu_num_2(3,:),'g*','linewidth',1);
text(0.2,100,num2str(cc_1_max(3)),'color',[178 34 34]./255,'FontSize',fs,'FontName','time','interpreter','latex') 
text(0.7,100,num2str(cc_2_max(3)),'color',[61 145 64]./255,'FontSize',fs,'FontName','time','interpreter','latex') 
xlim([0 1]);  ylim([-10 130])
xticks(0:0.2:1); yticks([0 100]);
xlabel('$Time \; (s)$','interpreter','latex'); ylabel('$\Omega_3 \; (rad \cdot s^{-1})$','FontSize',fs,'FontName','time','interpreter','latex')




subplot(3,2,2)
hold on; box on;
ax = gca; ax.TickLabelInterpreter = 'latex'; ax.FontSize = fs;
plot(t_exp,wu_exp(4,:),'k','linewidth',4);
plot(t_exp,wu_num_1(4,:),'r','linewidth',2);
plot(t_exp,wu_num_2(4,:),'g*','linewidth',1);
text(0.2,220,num2str(cc_1_max(4)),'color',[178 34 34]./255,'FontSize',fs,'FontName','time','interpreter','latex') 
text(0.7,220,num2str(cc_2_max(4)),'color',[61 145 64]./255,'FontSize',fs,'FontName','time','interpreter','latex') 
xlim([0 1]);  ylim([-140 300])
xticks(0:0.2:1); yticks([0 250]);
xticklabels([]); ylabel('$u_1 \; (\mu m \cdot s^{-1})$','FontSize',fs,'FontName','time','interpreter','latex')



subplot(3,2,4)
hold on; box on;
ax = gca; ax.TickLabelInterpreter = 'latex'; ax.FontSize = fs;
plot(t_exp,wu_exp(5,:),'k','linewidth',4);
plot(t_exp,wu_num_1(5,:),'r','linewidth',2);
plot(t_exp,wu_num_2(5,:),'g*','linewidth',1);
text(0.2,250,num2str(cc_1_max(5)),'color',[178 34 34]./255,'FontSize',fs,'FontName','time','interpreter','latex') 
text(0.7,250,num2str(cc_2_max(5)),'color',[61 145 64]./255,'FontSize',fs,'FontName','time','interpreter','latex') 
xlim([0 1]);  ylim([-350 350])
xticks(0:0.2:1); yticks([-300 0 300]);
xticklabels([]); ylabel('$u_2 \; (\mu m \cdot s^{-1})$','FontSize',fs,'FontName','time','interpreter','latex')


subplot(3,2,6)
hold on; box on;
ax = gca; ax.TickLabelInterpreter = 'latex'; ax.FontSize = fs;
plot(t_exp,wu_exp(6,:),'k','linewidth',4);
plot(t_exp,wu_num_1(6,:),'r','linewidth',2);
plot(t_exp,wu_num_2(6,:),'g*','linewidth',1);
text(0.2,200,num2str(cc_1_max(6)),'color',[178 34 34]./255,'FontSize',fs,'FontName','time','interpreter','latex') 
text(0.7,200,num2str(cc_2_max(6)),'color',[61 145 64]./255,'FontSize',fs,'FontName','time','interpreter','latex') 
xlim([0 1]);  ylim([-150 300])
xticks(0:0.2:1); yticks([0 250]);
xlabel('$Time \; (s)$','interpreter','latex'); ylabel('$u_3 \; (\mu m \cdot s^{-1})$','FontSize',fs,'FontName','time','interpreter','latex')

return

end


%% Plot the FFT power spectra.

if 0
fs=34;

 
figure(2)
clf; 
set(gcf,'color','w','units','normalized','outerposition',[0 0 1 1]) 

% Angular speeds
subplot(3,2,1)
hold on; box on;
ax = gca; ax.TickLabelInterpreter = 'latex'; ax.FontSize = fs;
plot(FD,PS_wu_exp(1,:),'k','linewidth',2)
scatter(Fre_wu_exp(1),PS_wu_exp(1,find(FD==Fre_wu_exp(1))),166,'r','filled')
text(Fre_wu_exp(1),PS_wu_exp(1,find(FD==Fre_wu_exp(1))),[num2str(Fre_wu_exp(1)),'HZ'],'fontname','times','fontsize',fs)
ax = gca; ax.TickLabelInterpreter = 'latex'; 
xlim([0 45]);
xticks([0 10 20 30 40]);
xticklabels([]); ylabel('$S_{\Omega_1}$','fontname','times','fontsize',fs+10,'interpreter','latex')



subplot(3,2,3)
hold on; box on;
ax = gca; ax.TickLabelInterpreter = 'latex'; ax.FontSize = fs;
plot(FD,PS_wu_exp(2,:),'k','linewidth',2)
scatter(Fre_wu_exp(2),PS_wu_exp(2,find(FD==Fre_wu_exp(2))),166,'r','filled')
text(Fre_wu_exp(2),PS_wu_exp(2,find(FD==Fre_wu_exp(2))),[num2str(Fre_wu_exp(2)),'HZ'],'fontname','times','fontsize',fs)
ax = gca; ax.TickLabelInterpreter = 'latex';
xlim([0 45]);
xticks([0 10 20 30 40]);
xticklabels([]); ylabel('$S_{\Omega_2}$','fontname','times','fontsize',fs+10,'interpreter','latex')


subplot(3,2,5)
hold on; box on;
ax = gca; ax.TickLabelInterpreter = 'latex'; ax.FontSize = fs;
plot(FD,PS_wu_exp(3,:),'k','linewidth',2)
scatter(Fre_wu_exp(3),PS_wu_exp(3,find(FD==Fre_wu_exp(3))),166,'r','filled')
text(Fre_wu_exp(3),PS_wu_exp(3,find(FD==Fre_wu_exp(3))),[num2str(Fre_wu_exp(3)),'HZ'],'fontname','times','fontsize',fs)
ax = gca; ax.TickLabelInterpreter = 'latex'; 
xlim([0 45]);
xticks([0 10 20 30 40]);
xlabel('Frequency domain','fontname','times','fontsize',fs+10); 
ylabel('$S_{\Omega_3}$','fontname','times','fontsize',fs+10,'interpreter','latex')


subplot(3,2,2)
hold on; box on;
ax = gca; ax.TickLabelInterpreter = 'latex'; ax.FontSize = fs;
plot(FD,PS_wu_exp(4,:),'k','linewidth',2)
scatter(Fre_wu_exp(4),PS_wu_exp(4,find(FD==Fre_wu_exp(4))),166,'r','filled')
text(Fre_wu_exp(4),PS_wu_exp(4,find(FD==Fre_wu_exp(4))),[num2str(Fre_wu_exp(4)),'HZ'],'fontname','times','fontsize',fs)
ax = gca; ax.TickLabelInterpreter = 'latex'; 
xlim([0 45]);
xticks([0 10 20 30 40]);
xlabel([]); ylabel('$S_{u_1}$','fontname','times','fontsize',fs+10,'interpreter','latex')


subplot(3,2,4)
hold on; box on;
ax = gca; ax.TickLabelInterpreter = 'latex'; ax.FontSize = fs;
plot(FD,PS_wu_exp(5,:),'k','linewidth',2)
scatter(Fre_wu_exp(5),PS_wu_exp(5,find(FD==Fre_wu_exp(5))),166,'r','filled')
text(Fre_wu_exp(5),PS_wu_exp(5,find(FD==Fre_wu_exp(5))),[num2str(Fre_wu_exp(5)),'HZ'],'fontname','times','fontsize',fs)
ax = gca; ax.TickLabelInterpreter = 'latex'; %ax.FontSize = 10; 
xlim([0 45]);
xticks([0 10 20 30 40]);
xlabel([]); ylabel('$S_{u_2}$','fontname','times','fontsize',fs+10,'interpreter','latex')


subplot(3,2,6)
hold on; box on;
ax = gca; ax.TickLabelInterpreter = 'latex'; ax.FontSize = fs;
plot(FD,PS_wu_exp(6,:),'k','linewidth',2)
scatter(Fre_wu_exp(6),PS_wu_exp(6,find(FD==Fre_wu_exp(6))),166,'r','filled')
text(Fre_wu_exp(6),PS_wu_exp(6,find(FD==Fre_wu_exp(6))),[num2str(Fre_wu_exp(6)),'HZ'],'fontname','times','fontsize',fs)
ax = gca; ax.TickLabelInterpreter = 'latex'; %ax.FontSize = 10; 
xlim([0 45]);
xticks([0 10 20 30 40]);
xlabel('Frequency domain','fontname','times','fontsize',fs+10); 
ylabel('$S_{u_3}$','fontname','times','fontsize',fs+10,'interpreter','latex')

return

end



%% Plot cross correlation.

if 0
fs=34;

figure(3)
clf; 
set(gcf,'color','w','units','normalized','outerposition',[0 0 1 1]) 
hold on; box on;
ax = gca; ax.TickLabelInterpreter = 'latex'; ax.FontSize = fs;
plot(lags_1(1,:),cc_1(1,:),'k','linewidth',2)
plot(lags_1(2,:),cc_1(2,:),'r','linewidth',2)
plot(lags_1(3,:),cc_1(3,:),'g','linewidth',2)
plot(lags_1(4,:),cc_1(4,:),'b','linewidth',2)
plot(lags_1(5,:),cc_1(5,:),'m','linewidth',2)
plot(lags_1(6,:),cc_1(6,:),'c','linewidth',2)

return

end

%% Plot sperm snapshot.


load 03_Head_Fixed_Frame_Method_Rotations_Traces_raw_5_Microns_Difraction_head_spin_rxmHeadNBT.mat;
b1 = head_tangent{sp};  % normalized. 3*nt.
b2 = head_normal{sp}; 
b3 = head_binormal{sp}; 


load 01_lab_frame_raw_traces_2017_2018_rxmSSppxyz.mat;
Z{sp}=Z{sp}-65;
arc_mat = arclength(X{sp},Y{sp},Z{sp});
arc_mean = mean(arc_mat(end,:));




%==============================================================================================================
% Generate sperm head: humanoid size.
    Ndiv = 60;
    radius_plot = 1.2/56; %normalized --> doesn't matter if head_type = 1.
    head_type = 1;    
    [Xhead,Yhead,Zhead] = sperm_head3D_shape(radius_plot, head_type, Ndiv);%size(Zhead)=[61 61]
    scal=1.2*arc_mean;  %In this plot, the head is for visualization, not used for reconstruction calculation, so the scale parameter is flexible.
    Xhead=Xhead*scal;
    Yhead=Yhead*scal;
    Zhead=Zhead*scal;
    Xneck=max(max(Xhead)); Xhead=Xhead-Xneck;
     
    M1=size(Xhead,1);M2=size(Xhead,2);M=M1*M2;
    % head upside
    [zup_row,zup_col]=find(Zhead>=0); %length(zup_row)=1891
    Zhead_up=nan(size(Zhead));Xhead_up=nan(size(Xhead));Yhead_up=nan(size(Yhead));
    for i = 1:length(zup_row)
        Zhead_up(zup_row(i),zup_col(i))=Zhead(zup_row(i),zup_col(i));
        Xhead_up(zup_row(i),zup_col(i))=Xhead(zup_row(i),zup_col(i));
        Yhead_up(zup_row(i),zup_col(i))=Yhead(zup_row(i),zup_col(i));
    end
    xup=[reshape(Xhead_up,M,1);reshape(Yhead_up,M,1);reshape(Zhead_up,M,1)]; 
    % head downside
    [zdown_row,zdown_col]=find(Zhead<=0);%length(zdown_row)=1891
    Zhead_down=nan(size(Zhead));Xhead_down=nan(size(Xhead));Yhead_down=nan(size(Yhead));
    for i = 1:length(zdown_row)
        Zhead_down(zdown_row(i),zdown_col(i))=Zhead(zdown_row(i),zdown_col(i));
        Xhead_down(zdown_row(i),zdown_col(i))=Xhead(zdown_row(i),zdown_col(i));
        Yhead_down(zdown_row(i),zdown_col(i))=Yhead(zdown_row(i),zdown_col(i));
    end   
    xdown=[reshape(Xhead_down,M,1);reshape(Yhead_down,M,1);reshape(Zhead_down,M,1)]; 

    

    %==============================================================================================================
    %Plot
    figure(4) 
    set(gcf,'color','w','units','normalized','outerposition',[0 0 1 1]) 
    clf;hold on;
    box off;axis equal;axis off;
    view(125,23)
    
    i_nt=10;
    
    % Plot tail. ===========================================================
    plot3(X{sp}(:,i_nt),Y{sp}(:,i_nt),Z{sp}(:,i_nt),'color',[36/255,100/255,171/255],'linewidth',4);
        
    % Plot head. ===========================================================
    B=[b1(:,i_nt),b2(:,i_nt),b3(:,i_nt)]; 
    neck=[X{sp}(1,i_nt),Y{sp}(1,i_nt),Z{sp}(1,i_nt)];            
    % head surface rotation and translation
    % upside head
    xupR=ApplyRotationMatrix(B,xup);
    xupT=TranslatePoints(xupR,neck);
    [head_up_x,head_up_y,head_up_z]=ExtractComponents(xupT);
    head_up_x=reshape(head_up_x,M1,M2);head_up_y=reshape(head_up_y,M1,M2);head_up_z=reshape(head_up_z,M1,M2);
    % downside head
    xdownR=ApplyRotationMatrix(B,xdown);
    xdownT=TranslatePoints(xdownR,neck);
    [head_down_x,head_down_y,head_down_z]=ExtractComponents(xdownT);
    head_down_x=reshape(head_down_x,M1,M2);head_down_y=reshape(head_down_y,M1,M2);head_down_z=reshape(head_down_z,M1,M2);
    % plot sperm head
    source_ligth = [90 80];
    k = [0.9 1 1 5]; 
    s1= surfl(head_up_x,head_up_y,head_up_z,source_ligth,k );
    s2= surfl(head_down_x,head_down_y,head_down_z,source_ligth,k );
    s1.EdgeColor='none' ;
    s1.FaceColor=  [0.6 0 0];
    s2.EdgeColor='none' ;
    s2.FaceColor=  [0 0.6 0];
  
    % Plot head basis vectors. ===========================================================
    head_a1 =  0.5*4.5; % semi-axis length of the longitudinal head axis.
    hc = [neck(1),neck(2),neck(3)] - b1(:,i_nt)'*head_a1; % head center coordinate.
    % Vector: line segment.
    eff_T=-9; eff_N=-9; eff_B=7; 
    quiver3(hc(1),hc(2),hc(3),eff_T*B(1,1),eff_T*B(2,1),eff_T*B(3,1),'k','linewidth',2); % arrow indicating head tangential vector
    quiver3(hc(1),hc(2),hc(3),eff_N*B(1,2),eff_N*B(2,2),eff_N*B(3,2),'k','linewidth',2); % arrow indicating head normal vector
    quiver3(hc(1),hc(2),hc(3),eff_B*B(1,3),eff_B*B(2,3),eff_B*B(3,3),'k','linewidth',2); % arrow indicating head binormal vector
    % Vector: arrowhead.
    eff_T_=-8; eff_N_=-8; eff_B_=6;
    [cone_T_x, cone_T_y, cone_T_z]=meshgrid(0:50:100, 0:50:100, -30:30:30);
    cone_T_u=repmat(-B(1,1),3,3,3);   
    cone_T_v=repmat(-B(2,1),3,3,3); 
    cone_T_w=repmat(-B(3,1),3,3,3); 
    cone_T_sxyz=[hc(1)+eff_T_*B(1,1),hc(2)+eff_T_*B(2,1),hc(3)+eff_T_*B(3,1)];
    hcone_T = coneplot(cone_T_x, cone_T_y, cone_T_z, cone_T_u, cone_T_v, cone_T_w,cone_T_sxyz(1),cone_T_sxyz(2),cone_T_sxyz(3),0.05);
    hcone_T.FaceColor = 'k';
    hcone_T.EdgeColor = 'none';
    hcone_T.DiffuseStrength = 0.8;  
    cone_N_u=repmat(-B(1,2),3,3,3);   
    cone_N_v=repmat(-B(2,2),3,3,3); 
    cone_N_w=repmat(-B(3,2),3,3,3); 
    cone_N_sxyz=[hc(1)+eff_N_*B(1,2),hc(2)+eff_N_*B(2,2),hc(3)+eff_N_*B(3,2)];
    hcone_N = coneplot(cone_T_x, cone_T_y, cone_T_z, cone_N_u, cone_N_v, cone_N_w,cone_N_sxyz(1),cone_N_sxyz(2),cone_N_sxyz(3),0.05);
    hcone_N.FaceColor = 'k';
    hcone_N.EdgeColor = 'none';
    hcone_N.DiffuseStrength = 0.8;  
    cone_B_u=repmat(B(1,3),3,3,3);   
    cone_B_v=repmat(B(2,3),3,3,3); 
    cone_B_w=repmat(B(3,3),3,3,3); 
    cone_B_sxyz=[hc(1)+eff_B_*B(1,3),hc(2)+eff_B_*B(2,3),hc(3)+eff_B_*B(3,3)];
    hcone_B = coneplot(cone_T_x, cone_T_y, cone_T_z, cone_B_u, cone_B_v, cone_B_w,cone_B_sxyz(1),cone_B_sxyz(2),cone_B_sxyz(3),0.05);
    hcone_B.FaceColor = 'k';
    hcone_B.EdgeColor = 'none';
    hcone_B.DiffuseStrength = 0.8;  
    
    
    
    % Plot angular speeds' concepts. ===========================================================
    [cone_ang_x, cone_ang_y, cone_ang_z]=meshgrid(-300:300:300, -300:300:300, -300:300:300);
    theta=0:pi/100:3/2*pi;
    % angular velocity around the tangential vector.
    cirT_y=3*sin(theta);
    cirT_z=3*cos(theta);
    cirT_x=cirT_y.*0;
    cone_angT_v=repmat(0,3,3,3);   
    cone_angT_w=repmat(1,3,3,3); 
    cone_angT_u=repmat(0,3,3,3); 
    cone_angT_sxyz=[cirT_x(end),cirT_y(end),cirT_z(end)];
    position_angT = [hc(1) hc(2) hc(3) ] + 1.1*[eff_T*B(1,1) eff_T*B(2,1) eff_T*B(3,1)];
    cirT=ApplyRotationMatrix(B,[cirT_x(:); cirT_y(:); cirT_z(:)]);
    cirT=TranslatePoints(cirT,position_angT);
    [cirT_x,cirT_y,cirT_z]=ExtractComponents(cirT);
    cone_angT=ApplyRotationMatrix(B,[cone_angT_u(:); cone_angT_v(:); cone_angT_w(:)]);
    [cone_angT_u, cone_angT_v, cone_angT_w]=ExtractComponents(cone_angT);
    cone_angT_u=reshape(cone_angT_u,3,3,3);  cone_angT_v=reshape(cone_angT_v,3,3,3);  cone_angT_w=reshape(cone_angT_w,3,3,3);
    cone_angT_sxyz=ApplyRotationMatrix(B,cone_angT_sxyz(:));
    cone_angT_sxyz=TranslatePoints(cone_angT_sxyz,position_angT);
    plot3(cirT_x,cirT_y,cirT_z,'m','linewidth',2);
    hcone_angT = coneplot(cone_ang_x, cone_ang_y, cone_ang_z, cone_angT_u, cone_angT_v, cone_angT_w,cone_angT_sxyz(1),cone_angT_sxyz(2),cone_angT_sxyz(3),0.01);
    hcone_angT.FaceColor = 'm';
    hcone_angT.EdgeColor = 'none';
    hcone_angT.DiffuseStrength = 0.8;  
    % angular velocity around the normal vector.
    cirN_x=3*sin(theta);
    cirN_z=3*cos(theta);
    cirN_y=cirN_x.*0;
    cone_angN_u=repmat(-1,3,3,3);   
    cone_angN_w=repmat(0,3,3,3); 
    cone_angN_v=repmat(0,3,3,3); 
    cone_angN_sxyz=[cirN_x(1),cirN_y(1),cirN_z(1)];
    position_angN = [hc(1) hc(2) hc(3)] + 1.1*[eff_N*B(1,2) eff_N*B(2,2) eff_N*B(3,2)];
    cirN=ApplyRotationMatrix(B,[cirN_x(:); cirN_y(:); cirN_z(:)]);
    cirN=TranslatePoints(cirN,position_angN);
    [cirN_x,cirN_y,cirN_z]=ExtractComponents(cirN);
    cone_angN=ApplyRotationMatrix(B,[cone_angN_u(:); cone_angN_v(:); cone_angN_w(:)]);
    [cone_angN_u, cone_angN_v, cone_angN_w]=ExtractComponents(cone_angN);
    cone_angN_u=reshape(cone_angN_u,3,3,3);  cone_angN_v=reshape(cone_angN_v,3,3,3);  cone_angN_w=reshape(cone_angN_w,3,3,3);
    cone_angN_sxyz=ApplyRotationMatrix(B,cone_angN_sxyz(:));
    cone_angN_sxyz=TranslatePoints(cone_angN_sxyz,position_angN);
    plot3(cirN_x,cirN_y,cirN_z,'m','linewidth',2);
    hcone_angN = coneplot(cone_ang_x, cone_ang_y, cone_ang_z, cone_angN_u, cone_angN_v, cone_angN_w,cone_angN_sxyz(1),cone_angN_sxyz(2),cone_angN_sxyz(3),0.01);
    hcone_angN.FaceColor = 'm';
    hcone_angN.EdgeColor = 'none';
    hcone_angN.DiffuseStrength = 0.8;  
    % angular velocity around the binormal vector.
    cirB_x=3*sin(theta);
    cirB_y=3*cos(theta);
    cirB_z=cirB_x.*0;
    cone_angB_u=repmat(-1,3,3,3);   
    cone_angB_v=repmat(0,3,3,3); 
    cone_angB_w=repmat(0,3,3,3); 
    cone_angB_sxyz=[cirB_x(1),cirB_y(1),cirB_z(1)];
    position_angB = [hc(1) hc(2) hc(3)] + 1.1*[eff_B*B(1,3) eff_B*B(2,3) eff_B*B(3,3)];
    cirB=ApplyRotationMatrix(B,[cirB_x(:); cirB_y(:); cirB_z(:)]);
    cirB=TranslatePoints(cirB,position_angB);
    [cirB_x,cirB_y,cirB_z]=ExtractComponents(cirB);
    cone_angB=ApplyRotationMatrix(B,[cone_angB_u(:); cone_angB_v(:); cone_angB_w(:)]);
    [cone_angB_u, cone_angB_v, cone_angB_w]=ExtractComponents(cone_angB);
    cone_angB_u=reshape(cone_angB_u,3,3,3);  cone_angB_v=reshape(cone_angB_v,3,3,3);  cone_angB_w=reshape(cone_angB_w,3,3,3);
    cone_angB_sxyz=ApplyRotationMatrix(B,cone_angB_sxyz(:));
    cone_angB_sxyz=TranslatePoints(cone_angB_sxyz,position_angB);
    plot3(cirB_x,cirB_y,cirB_z,'m','linewidth',2);
    hcone_angB = coneplot(cone_ang_x, cone_ang_y, cone_ang_z, cone_angB_u, cone_angB_v, cone_angB_w,cone_angB_sxyz(1),cone_angB_sxyz(2),cone_angB_sxyz(3),0.01);
    hcone_angB.FaceColor = 'm';
    hcone_angB.EdgeColor = 'none';
    hcone_angB.DiffuseStrength = 0.8;  


    
    % Plot linear speeds' concepts. ===========================================================
    % Vector: line segment.
    Ts = [hc(1),hc(2),hc(3)] + 1.2*[eff_T*B(1,1),eff_T*B(2,1),eff_T*B(3,1)]; % start from the tangential vector's distal end.
    Ns = [hc(1),hc(2),hc(3)] + 1.2*[eff_N*B(1,2),eff_N*B(2,2),eff_N*B(3,2)]; % start from the normal vector.
    Bs = [hc(1),hc(2),hc(3)] + 1.3*[eff_B*B(1,3),eff_B*B(2,3),eff_B*B(3,3)]; % start from the head binormal vector.
    eff_Ts=-4; eff_Ns=-4; eff_Bs=2; 
    quiver3(Ts(1),Ts(2),Ts(3),eff_Ts*B(1,1),eff_Ts*B(2,1),eff_Ts*B(3,1),'c','linewidth',2); % arrow indicating head tangential vector
    quiver3(Ns(1),Ns(2),Ns(3),eff_Ns*B(1,2),eff_Ns*B(2,2),eff_Ns*B(3,2),'c','linewidth',2); % arrow indicating head normal vector
    quiver3(Bs(1),Bs(2),Bs(3),eff_Bs*B(1,3),eff_Bs*B(2,3),eff_Bs*B(3,3),'c','linewidth',2); % arrow indicating head binormal vector
    % Vector: arrowhead.
    eff_Ts_=-4; eff_Ns_=-4; eff_Bs_=2;
    [cone_Ts_x, cone_Ts_y, cone_Ts_z]=meshgrid(0:50:100, 0:50:100, -30:30:30);
    cone_Ts_u=repmat(-B(1,1),3,3,3);   
    cone_Ts_v=repmat(-B(2,1),3,3,3); 
    cone_Ts_w=repmat(-B(3,1),3,3,3); 
    cone_Ts_sxyz=[Ts(1) Ts(2) Ts(3)] + [eff_Ts_*B(1,1) eff_Ts_*B(2,1) eff_Ts_*B(3,1)];
    hcone_Ts = coneplot(cone_Ts_x, cone_Ts_y, cone_Ts_z, cone_Ts_u, cone_Ts_v, cone_Ts_w,cone_Ts_sxyz(1),cone_Ts_sxyz(2),cone_Ts_sxyz(3),0.05);
    hcone_Ts.FaceColor = 'c';
    hcone_Ts.EdgeColor = 'none';
    hcone_Ts.DiffuseStrength = 0.8;  
    cone_Ns_u=repmat(-B(1,2),3,3,3);   
    cone_Ns_v=repmat(-B(2,2),3,3,3); 
    cone_Ns_w=repmat(-B(3,2),3,3,3); 
    cone_Ns_sxyz=[Ns(1)+eff_Ns_*B(1,2),Ns(2)+eff_Ns_*B(2,2),Ns(3)+eff_Ns_*B(3,2)];
    hcone_Ns = coneplot(cone_Ts_x, cone_Ts_y, cone_Ts_z, cone_Ns_u, cone_Ns_v, cone_Ns_w,cone_Ns_sxyz(1),cone_Ns_sxyz(2),cone_Ns_sxyz(3),0.05);
    hcone_Ns.FaceColor = 'c';
    hcone_Ns.EdgeColor = 'none';
    hcone_Ns.DiffuseStrength = 0.8;  
    cone_Bs_u=repmat(B(1,3),3,3,3);   
    cone_Bs_v=repmat(B(2,3),3,3,3); 
    cone_Bs_w=repmat(B(3,3),3,3,3); 
    cone_Bs_sxyz=[Bs(1)+eff_Bs_*B(1,3),Bs(2)+eff_Bs_*B(2,3),Bs(3)+eff_Bs_*B(3,3)];
    hcone_Bs = coneplot(cone_Ts_x, cone_Ts_y, cone_Ts_z, cone_Bs_u, cone_Bs_v, cone_Bs_w,cone_Bs_sxyz(1),cone_Bs_sxyz(2),cone_Bs_sxyz(3),0.05);
    hcone_Bs.FaceColor = 'c';
    hcone_Bs.EdgeColor = 'none';
    hcone_Bs.DiffuseStrength = 0.8;  
    
    
    % Plot Cartesian axes. ===========================================================
    O = [20 60 -8]; %original point.
    Cx=[-1; 0; 0]; Cy=[0; -1; 0]; Cz=[0; 0; 1]; C=[Cx Cy Cz];
    % Vector: line segment.
    eff_x=-4; eff_y=-4; eff_z=3; 
    quiver3(O(1),O(2),O(3),eff_x*C(1,1),eff_x*C(2,1),eff_x*C(3,1),'k','linewidth',2); % arrow indicating X axis.
    quiver3(O(1),O(2),O(3),eff_y*C(1,2),eff_y*C(2,2),eff_y*C(3,2),'k','linewidth',2); % arrow indicating Y axis.
    quiver3(O(1),O(2),O(3),eff_z*C(1,3),eff_z*C(2,3),eff_z*C(3,3),'k','linewidth',2); % arrow indicating Z axis.
    % Vector: arrowhead.
    eff_x_=-4; eff_y_=-4; eff_z_=3;
    [cone_X_x, cone_X_y, cone_X_z]=meshgrid(0:50:100, 0:50:100, -30:30:30);
    cone_X_u=repmat(-C(1,1),3,3,3);   
    cone_X_v=repmat(-C(2,1),3,3,3); 
    cone_X_w=repmat(-C(3,1),3,3,3); 
    cone_X_sxyz=[O(1)+eff_x_*C(1,1),O(2)+eff_x_*C(2,1),O(3)+eff_x_*C(3,1)];
    hcone_X = coneplot(cone_X_x, cone_X_y, cone_X_z, cone_X_u, cone_X_v, cone_X_w,cone_X_sxyz(1),cone_X_sxyz(2),cone_X_sxyz(3),0.05);
    hcone_X.FaceColor = 'k';
    hcone_X.EdgeColor = 'none';
    hcone_X.DiffuseStrength = 0.8;  
    cone_Y_u=repmat(-C(1,2),3,3,3);   
    cone_Y_v=repmat(-C(2,2),3,3,3); 
    cone_Y_w=repmat(-C(3,2),3,3,3); 
    cone_Y_sxyz=[O(1)+eff_y_*C(1,2),O(2)+eff_y_*C(2,2),O(3)+eff_y_*C(3,2)];
    hcone_Y = coneplot(cone_X_x, cone_X_y, cone_X_z, cone_Y_u, cone_Y_v, cone_Y_w,cone_Y_sxyz(1),cone_Y_sxyz(2),cone_Y_sxyz(3),0.05);
    hcone_Y.FaceColor = 'k';
    hcone_Y.EdgeColor = 'none';
    hcone_Y.DiffuseStrength = 0.8;  
    cone_Z_u=repmat(C(1,3),3,3,3);   
    cone_Z_v=repmat(C(2,3),3,3,3); 
    cone_Z_w=repmat(C(3,3),3,3,3); 
    cone_Z_sxyz=[O(1)+eff_z_*C(1,3),O(2)+eff_z_*C(2,3),O(3)+eff_z_*C(3,3)];
    hcone_Z = coneplot(cone_X_x, cone_X_y, cone_X_z, cone_Z_u, cone_Z_v, cone_Z_w,cone_Z_sxyz(1),cone_Z_sxyz(2),cone_Z_sxyz(3),0.05);
    hcone_Z.FaceColor = 'k';
    hcone_Z.EdgeColor = 'none';
    hcone_Z.DiffuseStrength = 0.8;
    
    
    
    camlight right
    lighting gouraud
    