function rxm_LabTrace_Plot

% This function plots the observed human sperm's trajectories in the lab frame, both the experimental one and representative reconstructed ones (with and without a boundary).

clc

sp=8;    

%% Plot the experimental trace.

if 0
    
    %================================================================================================================
    % Load and process the flagellar data.
    load 03_Head_Fixed_Frame_Method_Rotations_Traces_raw_5_Microns_Difraction_head_spin_rxmHeadNBT.mat;
    b1=head_tangent{sp};  % size(b1)=[3,nt]
    b2=head_normal{sp}; 
    b3=head_binormal{sp}; 
    clear X Y Z
    
    load 01_lab_frame_raw_traces_2017_2018_rxmSSppxyz.mat;
    Z{sp}=Z{sp}-65;
    ns = size(Z{sp},1);
    nt = size(Z{sp},2);
    for i_nt=1:nt
        arc_mat(:,i_nt) = arclength(X{sp}(:,i_nt), Y{sp}(:,i_nt), Z{sp}(:,i_nt));
    end
    arc_mean = mean(arc_mat(end,:)) % measured as 46.32.

    
    %================================================================================================================
    % Generate sperm head: humanoid size.
    Ndiv = 60;
    radius_plot = 1.2/56;  % normalized --> doesn't matter if head_type = 1.
    head_type = 1;    
    [Xhead,Yhead,Zhead] = sperm_head3D_shape(radius_plot, head_type, Ndiv); %size(Zhead)=[61 61]
    scal=1.2*arc_mean;  %In this sub-plot, the head is for visualization, not used for reconstruction calculation, so the scale parameter is flexible.
    Xhead=Xhead*scal;
    Yhead=Yhead*scal;
    Zhead=Zhead*scal;
    xneck=max(max(Xhead)); Xhead=Xhead-xneck;
    [~,np_head1] = max(Xhead(:));
    [~,np_head2] = min(Xhead(:));
    
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


    
    %================================================================================================================
    % Plot.
    figure(1)
    clf;
    set(gcf,'color','w','units','normalized'); 
    hold on;axis equal; box off;axis off;
    view(-30,15);

    
    % Plot sperm tail.
    clrs = colormap(winter(nt));
    for i_nt = 1:nt
        pT = plot3(X{sp}(:,i_nt),Y{sp}(:,i_nt),Z{sp}(:,i_nt),'color',clrs(i_nt,:),'linewidth',3) ;
        % transparency
        pT.Color(4)=0.07;
    end
    % Plot tail's mid-point trajectory.
    plot3(X{sp}(round(ns/2),:),Y{sp}(round(ns/2),:),Z{sp}(round(ns/2),:)-10,'color',[0.8500    0.3250    0.0980],'linewidth',1.8)    
    scatter3(X{sp}(round(ns/2),1),Y{sp}(round(ns/2),1),Z{sp}(round(ns/2),1)-10,140,'k','filled')    

    
    % Plot head
    load FreeSperm_Frequency.mat;
    % (time in seconds)  
    dt = 1/90; % Each frame is taken every (1/293) seconds => 1.00 s corresponds to 293 video data frames (framerate was 293fps).
    tRange = 0:dt:(nt-1)*dt;
    ncycles = tRange(end)*HeadSpin_Freq(sp);%number of revolution cycles in the lab frame, during the recording time
    dt_spin = round(nt/ncycles); %nspin=13; nt=272; ncycle_spin=20;
    for i_nt=[1,3*dt_spin+3,6*dt_spin+4,10*dt_spin+6,13*dt_spin+8,16*dt_spin+9,19*dt_spin+25]  % for sp#8
         B=[b1(:,i_nt),b2(:,i_nt),b3(:,i_nt)]; %orientation matrix (normalized) for sperm head
         neck=[X{sp}(2,i_nt),Y{sp}(2,i_nt),Z{sp}(2,i_nt)];     %neck=[X{sp}(1,i_nt),Y{sp}(1,i_nt),Z{sp}(1,i_nt)];            
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
         
         if i_nt==1
             head_down_z_initial = head_down_z;
         end
         
         plot3(X{sp}(:,i_nt),Y{sp}(:,i_nt),Z{sp}(:,i_nt),'color',clrs(i_nt,:),'linewidth',3) ;
    end
    

    % Plot boundary.
    % The boundary height is determined according to the initial position of head center.
    head_down_z_temp = head_down_z_initial(:);
    head_center_z = (head_down_z_temp(np_head1)+head_down_z_temp(np_head2))/2;
    % According to the experimental data, near-boundary cells have the boundary height about 5um.
    bound_height = 5;
    % set illustrated boundary size.
    XL=[min(X{sp}(:)) max(X{sp}(:))];
    YL=[min(Y{sp}(:)) max(Y{sp}(:))];
    ZL=[head_center_z-bound_height max(Z{sp}(:))];  
    %creating the side plane shading for 3D visualization.
    [xgrid,ygrid,zgrid] = meshgrid(linspace(XL(1),XL(2),10),...
                    linspace(YL(1),YL(2),10),...
                    linspace(ZL(1),ZL(2),10));                     
    CO(:,:,1) = ones(10)*.2;  %CO(:,:,1) = ones(10)*.0; % red
    CO(:,:,2) = ones(10)*.2; %CO(:,:,2) = ones(10)*.4470; % green
    CO(:,:,3) = ones(10)*0.2;  %CO(:,:,3) = ones(10)*0.7410; % blue
    % z-plane
    surfZ = surf(xgrid(:,:,1),ygrid(:,:,1),zgrid(:,:,1), CO,'FaceAlpha',0.1);  
    surfZ.EdgeColor='none';

    
    
    % Colorbar.
    colormap(winter); 
    swim_period = tRange(end); % The swimming time for sp8 is 3.0111s.
    swim_period = roundn(swim_period,-3);
    c = colorbar('southoutside','Ticks',[0,1],'TickLabels',{'0',num2str(swim_period)},'FontSize',28,'FontName','Times');  
    c.Label.String = 'Time (s)'; c.TickLabelInterpreter='latex'; 
    c.Position = [0.3 0.2 0.3 0.01];  % orientation: horizontal
    %c.Position = [0.85 0.2 0.01 0.7]; % orientation: vertical

    return
    
    if 1
        figure(3)
        clf;
        set(gcf,'color','w','units','normalized');
        hold on;axis equal;
        box off;axis off;
        view(54,-8);
        nrev = tRange(end)*CF_Freq(sp);
        dt_rev = round(nt/nrev);
        t_plot = 1:1*dt_rev+5;
        % Plot tail mid-point trajectory.
        plot3(X{sp}(round(ns/2),t_plot),Y{sp}(round(ns/2),t_plot),Z{sp}(round(ns/2),t_plot)-0.2,'color',[0.8500    0.3250    0.0980],'linewidth',10)    
        scatter3(X{sp}(round(ns/2),1),Y{sp}(round(ns/2),1),Z{sp}(round(ns/2),1)-0.2,2500,'k','filled')
    end
end



%% Plot the numerically reconstructed trace.

if 1
    
    %================================================================================================================
    % Load and process the flagellar data.
    load 01_lab_frame_raw_traces_2017_2018_rxmSSppxyz.mat; 
    nt_exp = size(Z{sp},2);
    dt=1/90;
    t_exp = 0:dt:(nt_exp-1)*dt;
    clearvars -except nt_exp t_exp sp

    
    %load XNodes_sp8_NoBoundStokeslet.mat;  %//////////////////////////////////////////////////////////////////
    load XNodes_sp8_BoundBlakelet_H0.5.mat;  %//////////////////////////////////////////////////////////////////
    X = Xs(Nhh+1:Nhh+Ns,:); 
    Y = Xs(2*Nhh+Ns+1:2*Nhh+2*Ns,:); 
    Z = Xs(3*Nhh+2*Ns+1:3*Nhh+3*Ns,:);
    ns=size(X,1);
    nt_num=size(X,2);
    head_tangent=head_tangent';
    head_normal=head_normal';
    head_binormal=nan(size(head_tangent)); % 3*nt
    for i_nt=1:nt_num
        head_binormal(:,i_nt)=cross(head_tangent(:,i_nt),head_normal(:,i_nt));          
    end
  
     
    % Time interpolation to get the reconstructed flagellum positions at the experimentally observed time points.
    s_inp = linspace(0,1,ns);
    t_dim = t_nd/(2*pi)/Freq_WF;
    [T_exp,S_exp] = meshgrid(t_exp,s_inp) ;
    [T_num,S_num] = meshgrid(t_dim,s_inp) ; 
    X = interp2(T_num,S_num,X,T_exp,S_exp,'spline'); 
    Y = interp2(T_num,S_num,Y,T_exp,S_exp,'spline'); 
    Z = interp2(T_num,S_num,Z,T_exp,S_exp,'spline');
    % Time interpolation to get the reconstructed head orientations at the experimentally observed time points.
    [T_exp,S_exp] = meshgrid(t_exp,[1 2 3]) ;
    [T_num,S_num] = meshgrid(t_dim,[1 2 3]) ; 
    head_tangent = interp2(T_num,S_num,head_tangent,T_exp,S_exp,'spline');             
    head_normal  = interp2(T_num,S_num,head_normal,T_exp,S_exp,'spline');             
    head_binormal = interp2(T_num,S_num,head_binormal,T_exp,S_exp,'spline');  

    
    % Spatial smooth of the flagellar shape, to remove the noise arising from the raw experimental data being reconstructed.
    ppx = 0.999; % 'ppx=1' means interpolation without cubic spline.
    ppy = 0.999;
    ppz = 0.999;  
    for i_nt = 1:length(t_exp)
        s_temp = arclength(X(:,i_nt), Y(:,i_nt), Z(:,i_nt)); 
        X(:,i_nt) = fnval(csaps(s_temp,X(:,i_nt),ppx),s_temp);
        Y(:,i_nt) = fnval(csaps(s_temp,Y(:,i_nt),ppy),s_temp);
        Z(:,i_nt) = fnval(csaps(s_temp,Z(:,i_nt),ppz),s_temp);
    end

    clearvars -except nt_exp ns X Y Z head_tangent head_normal head_binormal sp

    


    %================================================================================================================
    % Generate sperm head: simplified ellipsoid. 
    head_a1=2.0/45;
    head_a2=1.6/45;
    head_a3=1.0/45;
    nth=20;
    nphi=40;
    a=1;
    M=nth*nphi;

    %head: up side
    [xup1,xup2,xup3]=GenerateSphereSurfaceForVisualisation_up(nth,nphi,a);
    [M1 M2]=size(xup1);
    xup1=xup1*head_a1; xup1=xup1-max(xup1(:)); %这条代码之前没有添加，导致trace不完全正确，recon traces（包括正文和SI）都要重新画！
    xup2=xup2*head_a2;
    xup3=xup3*head_a3;
    xup=[reshape(xup1,M,1);reshape(xup2,M,1);reshape(xup3,M,1)]; 
    %head: down side
    [xdown1,xdown2,xdown3]=GenerateSphereSurfaceForVisualisation_down(nth,nphi,a);
    xdown1=xdown1*head_a1; xdown1=xdown1-max(xdown1(:)); %这条代码之前没有添加，导致trace不完全正确，recon traces（包括正文和SI）都要重新画！
    xdown2=xdown2*head_a2;
    xdown3=xdown3*head_a3;
    xdown=[reshape(xdown1,M,1);reshape(xdown2,M,1);reshape(xdown3,M,1)]; 
   
    [~,np_head1] = max(xdown1(:));
    [~,np_head2] = min(xdown1(:));

   
    
    %================================================================================================================
    % Plot 
    figure(2)
    clf;
    set(gcf,'color','w','units','normalized');
    hold on;axis equal;
    box off;axis off;
    view(-30,15);
    
    
    % Plot sperm tail
    clrs = colormap(winter(nt_exp));
    for i_nt = 1:nt_exp       
        pT = plot3(X(:,i_nt),Y(:,i_nt),Z(:,i_nt),'color',clrs(i_nt,:),'linewidth',3) ;
        % transparency.
        pT.Color(4)=0.07;            
    end    
    % Plot tail mid-point trajectory.
    plot3(X(round(ns/2),:),Y(round(ns/2),:),Z(round(ns/2),:)-0.2,'color',[0.8500    0.3250    0.0980],'linewidth',1.8)    
    scatter3(X(round(ns/2),1),Y(round(ns/2),1),Z(round(ns/2),1)-0.2,140,'k','filled')

    
    % plot sperm head
    load FreeSperm_Frequency.mat;
    % (time in seconds)  
    dt = 1/90; % Each frame is taken every (1/293) seconds => 1.00 s corresponds to 293 video data frames (framerate was 293fps).
    tRange = 0:dt:(nt_exp-1)*dt;
    ncycles = tRange(end)*HeadSpin_Freq(sp);%number of revolution cycles in the lab frame, during the recording time
    dt_spin = round(nt_exp/ncycles); %nspin=13; nt=272; ncycle_spin=20;
    for i_nt=[1,3*dt_spin+3,6*dt_spin+4,10*dt_spin+6,13*dt_spin+8,16*dt_spin+9,19*dt_spin+25]  % for sp#8
        B=[head_tangent(:,i_nt),head_normal(:,i_nt),head_binormal(:,i_nt)]; 
        neck=[X(2,i_nt),Y(2,i_nt),Z(2,i_nt)];          %neck=[X(1,i_nt),Y(1,i_nt),Z(1,i_nt)];            
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
        s1.FaceColor=  [0 0.6 0];
        s2.EdgeColor='none' ;
        s2.FaceColor=  [0.6 0 0];             
               
        if i_nt==1
            head_down_z_initial = head_down_z;
        end
        
        plot3(X(:,i_nt),Y(:,i_nt),Z(:,i_nt),'color',clrs(i_nt,:),'linewidth',3) ;
    end
    
    
    if 1%//////////////////////////////////////////////////////////////////
    % Plot boundary
    % find the boundary height basis: initial position of head center.
    head_down_z_temp = head_down_z_initial(:);
    head_center_z = (head_down_z_temp(np_head1)+head_down_z_temp(np_head2))/2;
    % The boundary height for different numerical cases varies.
    bound_height = 0.5;
    % set illustrated boundary size.
    XL=[min(X(:)) max(X(:))];
    YL=[min(Y(:)) max(Y(:))];
    ZL=[head_center_z-bound_height max(Z(:))];  
    %creating the side plane shading for 3D visualization======================
    [xgrid,ygrid,zgrid] = meshgrid(linspace(XL(1),XL(2),10),...
                    linspace(YL(1),YL(2),10),...
                    linspace(ZL(1),ZL(2),10));                     
    CO(:,:,1) = ones(10)*.2;  %CO(:,:,1) = ones(10)*.0; % red
    CO(:,:,2) = ones(10)*.2; %CO(:,:,2) = ones(10)*.4470; % green
    CO(:,:,3) = ones(10)*0.2;  %CO(:,:,3) = ones(10)*0.7410; % blue
    % z-plane
    surfZ = surf(xgrid(:,:,1),ygrid(:,:,1),zgrid(:,:,1), CO,'FaceAlpha',0.1);   %surf(xgrid(:,:,1),ygrid(:,:,1),zgrid(:,:,1)*0+head_center_z-bound_height, CO,'FaceAlpha',0.05); 
    surfZ.EdgeColor='none';
    end
    
    
    if 1
        figure(3)
        clf;
        set(gcf,'color','w','units','normalized');
        hold on;axis equal;
        box off;axis off;
        view(54,-8);
        nrev = tRange(end)*CF_Freq(sp);
        dt_rev = round(nt_exp/nrev);
        t_plot = 1:1*dt_rev+5;
        % Plot tail mid-point trajectory.
        plot3(X(round(ns/2),t_plot),Y(round(ns/2),t_plot),Z(round(ns/2),t_plot)-0.2,'color',[0.8500    0.3250    0.0980],'linewidth',10)    
        scatter3(X(round(ns/2),1),Y(round(ns/2),1),Z(round(ns/2),1)-0.2,2500,'k','filled')
    end
    
    return
end




%% Plot colorbar.

figure(4)
clf; set(gcf,'color','w','units','normalized');
hold on;axis equal; box off;axis off;
        
 
    colormap(winter); 
    swim_period = 3.011; %tRange(end); % The swimming time for sp8 is 3.0111s.
    swim_period = roundn(swim_period,-3);
    c = colorbar('eastoutside','Ticks',[0,1],'TickLabels',{'0',num2str(swim_period)},'FontSize',48,'FontName','Times');  
    c.Label.String = 'Time (s)'; c.TickLabelInterpreter='latex'; 
    %c.Position = [0.3 0.2 0.3 0.01];  % orientation: horizontal
    c.Position = [0.6 0.2 0.008 0.7]; % orientation: vertical

