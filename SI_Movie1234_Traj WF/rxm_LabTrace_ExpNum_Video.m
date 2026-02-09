function rxm_LabTrace_ExpNum_Video


% This function generates videos for the experimental sperm trace in lab frame, and
% compare it with the corresponding numerical reconstruction, with boundary considered.
% The boundary height is determined by the distance from the initial
% position of sperm head center.

clc

sp=8;

Video_exp=1;
Video_num=1;


%% Experimental trace. The boundary height is 5um, about 5/50=0.1L.

if Video_exp==1
    
    load 03_Head_Fixed_Frame_Method_Rotations_Traces_raw_5_Microns_Difraction_head_spin_rxmHeadNBT.mat;  
    % Note that the basis vectors (b1,b2,b3) here have already been normalized!
    b1=head_tangent{sp}; % size(b1)=[3,nt]
    b2=head_normal{sp}; % B=[b1,b2,b3]. The sequency of 'b1,2,3' is determined by head coordinates below (xup,xdown)
    b3=head_binormal{sp};% because 'xup=[head_a1;head_a2;head_a3]', the basis vectors should keep consistent, that is, 'b1' along 'head_a1', 'b2' along 'head_a2',...
    nt = size(b1,2);
    clear X Y Z
    
    
    load 01_lab_frame_raw_traces_2017_2018_rxmSSppxyz.mat;
    Z{sp}=Z{sp}-Z{sp}(1,1);
    %Z{sp}=Z{sp}-65;%This is to make the cell near the boundary look like swim close to the boundary in the video.
    ns = size(Z{sp},1);
    for i_nt=1:nt
        arc_mat(:,i_nt) = arclength(X{sp}(:,i_nt), Y{sp}(:,i_nt), Z{sp}(:,i_nt));
    end
    arc_mean = mean(arc_mat(end,:)); % measured as 46.32.


    
    %% generate sperm head: humanoid size

    %creating mesh for sperm head provide by Hermes
    Ndiv = 60;
    radius_plot = 1.2/56;%normalized
    head_type = 1;    
    [Xhead,Yhead,Zhead] = sperm_head3D_shape(radius_plot, head_type, Ndiv);%size(Zhead)=[61 61]
    scal=1.2*arc_mean;  %In this sub-plot, the head is for visualization, not used for reconstruction calculation, so the scale parameter is flexible.
    Xhead=Xhead*scal;Yhead=Yhead*scal;Zhead=Zhead*scal;
    % translate the head to make the 'neck' point at origin [0 0 0].
    xneck=max(max(Xhead));Xhead=Xhead-xneck;
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

    
    
    %% plot raw lab trajectory, 'ns' space points and 'nt' time instants
    % Font options for printing figures
    fs=32;
    fn='times';
    fig = 1;

    figure(fig) 

    v=VideoWriter(['LabTraj_Exp'],'MPEG-4');
    v.FrameRate=6;
    open(v);     

    % Locate the boundary, with the boundary distance to the initial head center position 5um.                      
    % Define head geometory.
    % neck (head centre) position to locate/ translate the sperm head
    B=[b1(:,1),b2(:,1),b3(:,1)]; % (initial) orientation matrix (normalized) for sperm head
    neck=[X{sp}(1,1),Y{sp}(1,1),Z{sp}(1,1)];          
    % head surface rotation and translation
    % downside head
    xdownR=ApplyRotationMatrix(B,xdown);
    xdownT=TranslatePoints(xdownR,neck);
    [~,~,head_down_z]=ExtractComponents(xdownT);
    head_down_z_initial=reshape(head_down_z,M1,M2);
    % Locate the boundary.                      
    head_down_z_temp = head_down_z_initial(:);
    head_center_z = (head_down_z_temp(np_head1)+head_down_z_temp(np_head2))/2;
    % According to the experimental data, near-boundary cells have boundary
    % height about 5um.
    bound_height = 5+2;  %"+2" because otherwise the head looks too close to the boundary surface.
    ZL_min = head_center_z-bound_height; 
        

    for i_nt= 11:nt-1 
        set(gcf,'color','w','units','normalized','outerposition',[0 0 1 1]) % make it white background
        clf;hold on;
        box on;axis equal;
        view(-35,20); % for sp#8   
        XL = [20 130];YL = [-5 70];ZL = [ZL_min 15];  
        
        % Plot the boundary.==========================================================               
        [xgrid,ygrid,zgrid] = meshgrid(linspace(XL(1),XL(2),10),...
                    linspace(YL(1),YL(2),10),...
                    linspace(ZL(1),ZL(2),10));                     
        CO_b(:,:,1) = ones(10)*1; 
        CO_b(:,:,2) = ones(10)*1; 
        CO_b(:,:,3) = ones(10)*1; 
        surf(xgrid(:,:,1),ygrid(:,:,1),zgrid(:,:,1)*0+ZL_min, CO_b),
        shading interp

        
        % plot sperm tail=========================================================
        % Transparency parameter.
        if i_nt>11
        for k=1:i_nt-11
            PlotTailCylinder_Movie_LabTrajProj(X{sp}(:,k),Y{sp}(:,k),Z{sp}(:,k),0.1,fig, XL(2),YL(2),ZL(1))
        end  
        end
        FAlpha_mat_temp = [-0.9:0.1:-0.5 -0.45:0.05:-0.25 0]; %FAlpha_mat_temp = linspace(-1,0,11);
        FAlpha_mat = 10.^FAlpha_mat_temp;               
        i_alpha=0;
        for k=i_nt-10:i_nt
            i_alpha = i_alpha+1;
            PlotTailCylinder_Movie_LabTraj(X{sp}(:,k),Y{sp}(:,k),Z{sp}(:,k),FAlpha_mat(i_alpha),fig, XL(2),YL(2),ZL(1))
        end
        %
        % plot tail mid-point trajectory and projection=========================================================
        plot3(X{sp}(round(ns/2),1:i_nt),Y{sp}(round(ns/2),1:i_nt),Z{sp}(round(ns/2),1:i_nt),'color',[0.8500    0.3250    0.0980],'linewidth',1.5)    
        plot3(X{sp}(round(ns/2),1:i_nt),Y{sp}(round(ns/2),1:i_nt),Z{sp}(round(ns/2),1:i_nt)*0+ZL(1)+0.1,'color',[0.8500    0.3250    0.0980],'linewidth',0.1) 
        scatter3(X{sp}(round(ns/2),i_nt),Y{sp}(round(ns/2),i_nt),Z{sp}(round(ns/2),i_nt),60,[0.8500    0.3250    0.0980],'filled')    
        scatter3(X{sp}(round(ns/2),i_nt),Y{sp}(round(ns/2),i_nt),Z{sp}(round(ns/2),i_nt)*0+ZL(1)+0.1,30,[0.8500    0.3250    0.0980],'filled') 
   
           
        
        % Plot head==================================================================
        % Define head geometory.
        % neck (head centre) position to locate/ translate the sperm head
        B=[b1(:,i_nt),b2(:,i_nt),b3(:,i_nt)]; %orientation matrix (normalized) for sperm head
        neck=[X{sp}(2,i_nt),Y{sp}(2,i_nt),Z{sp}(2,i_nt)];          
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
        % 3D head.
        source_ligth = [90 80];
        k = [0.9 1 1 5]; 
        s1= surfl(head_up_x,head_up_y,head_up_z,source_ligth,k );
        s2= surfl(head_down_x,head_down_y,head_down_z,source_ligth,k );
        s1.EdgeColor='none' ;
        s1.FaceColor=  [177 24 45]/255;  %[0.6 0 0];  
        s2.EdgeColor='none' ;
        s2.FaceColor=  [218 207 168]/255;  %[0 0.6 0];
        % head projection on XY plane.
        ss1= surf(head_up_x,head_up_y,head_up_z.*0+ZL(1) );
        ss2= surf(head_down_x,head_down_y,head_down_z.*0+ZL(1));
        cc = nan(size(s1.CData,1),size(s1.CData,2),3); cc(:,:,1)=0; cc(:,:,2)=0; cc(:,:,3)=0;ss1.CData=cc;
        cc = nan(size(s1.CData,1),size(s1.CData,2),3); cc(:,:,1)=0; cc(:,:,2)=0; cc(:,:,3)=0;ss2.CData=cc;             
        % head projection on XZ plane.
        ss1= surf(head_up_x,head_up_y.*0+YL(2),head_up_z);
        ss2= surf(head_down_x,head_down_y.*0+YL(2),head_down_z);
        cc = nan(size(s1.CData,1),size(s1.CData,2),3); cc(:,:,1)=0; cc(:,:,2)=0; cc(:,:,3)=0;ss1.CData=cc;
        cc = nan(size(s1.CData,1),size(s1.CData,2),3); cc(:,:,1)=0; cc(:,:,2)=0; cc(:,:,3)=0;ss2.CData=cc;             
        % head projection on YZ plane.
        ss1= surf(head_up_x.*0+XL(2),head_up_y,head_up_z);
        ss2= surf(head_down_x.*0+XL(2),head_down_y,head_down_z);
        cc = nan(size(s1.CData,1),size(s1.CData,2),3); cc(:,:,1)=0; cc(:,:,2)=0; cc(:,:,3)=0;ss1.CData=cc;
        cc = nan(size(s1.CData,1),size(s1.CData,2),3); cc(:,:,1)=0; cc(:,:,2)=0; cc(:,:,3)=0;ss2.CData=cc;             
             
        
        % Plot head basis==================================================================
        if 0
        % head center position (instantaneous)
        head_down_x_temp = head_down_x(:);
        head_down_y_temp = head_down_y(:);
        head_down_z_temp = head_down_z(:);
        head_center_x = (head_down_x_temp(np_head1)+head_down_x_temp(np_head2))/2;        
        head_center_y = (head_down_y_temp(np_head1)+head_down_y_temp(np_head2))/2;
        head_center_z = (head_down_z_temp(np_head1)+head_down_z_temp(np_head2))/2;
        O = [head_center_x, head_center_y, head_center_z]; 
        % head orientation (instantaneous)
        Cx=[-10; 0; 0]; 
        Cy=[0; -10; 0]; 
        Cz=[0; 0; 10]; 
        C_raw=[Cx; Cy; Cz];
        C_temp=ApplyRotationMatrix(B,C_raw);
        [C_new_x,C_new_y,C_new_z]=ExtractComponents(C_temp);
        C_new = [C_new_x,C_new_y,C_new_z];
        %
        Plot_Cartesian_Axes(fig,O,C_new') 
        % axes texts.
        text_x=[-12;0;0];
        text_y=[0;-12;0];
        text_z=[0;0;12];
        text_po=[text_x;text_y;text_z]; % position of text.
        text_R=ApplyRotationMatrix(B,text_po);
        text_RT=TranslatePoints(text_R,O);
        [text_RTx,text_RTy,text_RTz]=ExtractComponents(text_RT);
        text(text_RTx(1),text_RTy(1),text_RTz(1),'$e_3$','Fontname',fn,'FontSize',fs,'Interpreter','latex','fontweight','bold');
        text(text_RTx(2),text_RTy(2),text_RTz(2),'$e_1$','Fontname',fn,'FontSize',fs,'Interpreter','latex','fontweight','bold');
        text(text_RTx(3),text_RTy(3),text_RTz(3),'$e_2$','Fontname',fn,'FontSize',fs,'Interpreter','latex','fontweight','bold');
        end
        
        
        % Title: time.==================================================================      
        % (time in seconds)  
        dt = 1/90; % Each frame is taken every (1/90) seconds => 1.00 s corresponds to 90 video data frames (framerate was 90fps).
        inst = dt*(i_nt-1);
        format short
        inst = roundn(inst,-2);
        title(['Time: ',num2str(inst),'s'],'Fontname','Times New Roman','fontsize',fs,'fontweight','normal');
    
    
        xlim([XL(1),XL(2)]);ylim([YL(1),YL(2)]);zlim([ZL(1),ZL(2)]);
        xticks([0:20:140]);yticks([0:20:140]);zticks([-10:10:10])  
        ax = gca;ax.FontSize = fs; ax.TickLabelInterpreter = 'latex';
        xlh=xlabel('$x (\mu m)$ ','interpreter','latex','FontSize',1.5*fs);
        ylh=ylabel('$y (\mu m)$ ','interpreter','latex','FontSize',1.5*fs);
        zlh=zlabel('$z (\mu m)$ ','interpreter','latex','FontSize',1.5*fs,'Rotation',0);
        set(xlh,'position',[80 YL(1)-5 ZL(1)]-5);
        set(ylh,'position',[XL(1)-2 30 ZL(1)-5]);
        set(zlh,'position',[XL(1)-11 YL(2)+5 2]);%set(zlh,'position',[XL(1)-1 YL(2)+10 75]);
        
        
        l1=light;
        l1.Color = [1 1 1];
        l1.Style = 'local';  %'infinite';  %
        l1.Position = [30 60 20];
        l2=light;
        l2.Color = [1 1 1];
        l2.Style = 'local';  %'infinite';  %
        l2.Position = [80 35 20];
        l3=light;
        l3.Color = [1 1 1];
        l3.Style = 'local';  %'infinite';  %
        l3.Position = [120 5 20];
        lighting gouraud   %phong  %flat  %
        material metal  %shiny     


        drawnow;
        frame(i_nt)=getframe(gcf);
        writeVideo(v,frame(i_nt));
   
    end
    
    %return
end





%% Numerical reconstructed trace, with/ without a boundary.

if Video_num==1
    
    % Step 1: read numerical data (fine time discretization), and then get coarse-time interpolation and spatial smooth.
    load 01_lab_frame_raw_traces_2017_2018_rxmSSppxyz.mat; 
    nt = size(Z{sp},2);
    % (time in seconds)  
    dt=1/90; % Each frame is taken every (1/90) seconds => 1.00 s corresponds to 90 video data frames (framerate was 90fps).
    t_exp = 0:dt:(nt-1)*dt;
    for i_nt=1:nt
        arc_mat(:,i_nt) = arclength(X{sp}(:,i_nt), Y{sp}(:,i_nt), Z{sp}(:,i_nt));
    end
    arc_mean = mean(arc_mat(end,:)); % measured as 46.32.
    clearvars -except nt t_exp sp arc_mean 
    
    
    load XNodes_sp8_NoBoundStokeslet.mat;  %//////////////////////////////////////////////////////////////////
    %load XNodes_sp8_BoundBlakelet_H0.2.mat;    %//////////////////////////////////////////////////////////////////  
    X = Xs(Nhh+1:Nhh+Ns,:); Y = Xs(2*Nhh+Ns+1:2*Nhh+2*Ns,:); Z = Xs(3*Nhh+2*Ns+1:3*Nhh+3*Ns,:);
    nt_fine=size(X,2); %corresponding to the fine time discretization 'tRange' in the numerical reconstruction data, fine time for later flow field calcualtion.
    ns=size(X,1);
    head_tangent=head_tangent';
    head_normal=head_normal';
    head_binormal=nan(size(head_tangent)); % 3*nt
    for i_nt=1:nt_fine
        head_binormal(:,i_nt)=cross(head_tangent(:,i_nt),head_normal(:,i_nt));          
    end
    % Coarse-time interpolation.
    s_inp = linspace(0,1,ns);
    t_dim = t_nd/(2*pi)/Freq_WF;
    [Tcoarse,Scoarse] = meshgrid(t_exp,s_inp) ;% query points for coarser-time sperm tail coordinates.
    [Tfine,Sfine] = meshgrid(t_dim,s_inp) ; 
    X = interp2(Tfine,Sfine,X,Tcoarse,Scoarse,'spline'); % component 'x' of the tail coordinates, along arc length, at the instants of experimental recoarding 't_exp'.
    Y = interp2(Tfine,Sfine,Y,Tcoarse,Scoarse,'spline'); % component 'y' of the tail coordinates, along arc length, at the instants of experimental recoarding 't_exp'.
    Z = interp2(Tfine,Sfine,Z,Tcoarse,Scoarse,'spline'); % component 'z' of the tail coordinates, along arc length, at the instants of experimental recoarding 't_exp'.
    [Tcoarse,Scoarse] = meshgrid(t_exp,[1 2 3]) ;% query points for coarser-time sperm tail coordinates.
    [Tfine,Sfine] = meshgrid(t_dim,[1 2 3]) ; 
    head_tangent = interp2(Tfine,Sfine,head_tangent,Tcoarse,Scoarse,'spline');             
    head_normal  = interp2(Tfine,Sfine,head_normal,Tcoarse,Scoarse,'spline');             
    head_binormal = interp2(Tfine,Sfine,head_binormal,Tcoarse,Scoarse,'spline'); 
    b1 = head_tangent;
    b2 = head_normal;
    b3 = head_binormal;
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
    clearvars -except nt ns X Y Z b1 b2 b3 sp arc_mean ncycles 


    
    % Step 2: Generate sperm head: ellipsoid size semi-axes. 
    % The values can be modified, but their ratios should be unchanged.
    head_a1=1*arc_mean* 2.0/45; %
    head_a2=1*arc_mean* 1.6/45; %
    head_a3=1*arc_mean* 1.0/45; %
    %head:up+down£¬different colors,surf
    nth=20; 
    nphi=40; 
    a=1;
    M=nth*nphi;
    %head: up side
    [xup1,xup2,xup3]=GenerateSphereSurfaceForVisualisation_up(nth,nphi,a);
    [M1 M2]=size(xup1);
    xup1=xup1*head_a1; xup1=xup1-max(xup1(:));
    xup2=xup2*head_a2;
    xup3=xup3*head_a3;
    xup=[reshape(xup1,M,1);reshape(xup2,M,1);reshape(xup3,M,1)]; 
    %head: down side
    [xdown1,xdown2,xdown3]=GenerateSphereSurfaceForVisualisation_down(nth,nphi,a);
    xdown1=xdown1*head_a1; xdown1=xdown1-max(xdown1(:));
    xdown2=xdown2*head_a2;
    xdown3=xdown3*head_a3;
    xdown=[reshape(xdown1,M,1);reshape(xdown2,M,1);reshape(xdown3,M,1)]; 
    % locate the head axis for later location of head center in the lab frame.
    [~,np_head1] = max(xdown1(:));
    [~,np_head2] = min(xdown1(:));


    
        
    %% plot reconstructed lab trajectory, 'ns' space points and 'nt' time instants
    % Font options for printing figures
    fs=32;
    fn='times';
    fig = 2;
    
    figure(fig) 

    %v=VideoWriter(['LabTraj_NumBound02'],'MPEG-4');  %//////////////////////////////////////////////////////////////////
    v=VideoWriter(['LabTraj_NumNoBound'],'MPEG-4');  %//////////////////////////////////////////////////////////////////
    v.FrameRate=6;
    open(v);     

    % Locate the boundary, with the boundary distance to the initial head center position 0.2L/ infinite.                      
    % Define head geometory.
    % neck (head centre) position to locate/ translate the sperm head
    B=[b1(:,1),b2(:,1),b3(:,1)]; % (initial) orientation matrix (normalized) for sperm head
    neck=[X(1,1),Y(1,1),Z(1,1)];          
    % head surface rotation and translation
    % downside head
    xdownR=ApplyRotationMatrix(B,xdown);
    xdownT=TranslatePoints(xdownR,neck);
    [~,~,head_down_z]=ExtractComponents(xdownT);
    head_down_z_initial=reshape(head_down_z,M1,M2);
    % Locate the boundary.                      
    head_down_z_temp = head_down_z_initial(:);
    head_center_z = (head_down_z_temp(np_head1)+head_down_z_temp(np_head2))/2;
    % Set the boundary height.
    bound_height = 10;  %0.2+0.15;  %//////////////////////////////////////////////////////////////////
    % Note that for boundary height equal to 0.2L case, the visualized
    % boundary surface position is lower in Z coordinates, otherwise the
    % sperm body will contact with the boundary.
    ZL_min = head_center_z-bound_height; 
    
       
    X = X.*arc_mean; Y = Y.*arc_mean; Z = Z.*arc_mean;
    ZL_min = ZL_min*arc_mean;
    for i_nt= 11:nt-1  
        set(gcf,'color','w','units','normalized','outerposition',[0 0 1 1]) % make it white background
        clf;hold on;
        box on;axis equal;
        view(-35,20); % for sp#8          
        XL = [-0.9 1]*arc_mean;YL = [-0.8 0.5]*arc_mean;
        %ZL = [ZL_min 0.5*arc_mean]; % for 1 boundary %//////////////////////////////////////////////////////////////////
        ZL = [-0.3 0.25]*arc_mean; % for no boundary  %//////////////////////////////////////////////////////////////////
        
                       
        % Plot the boundary.==========================================================               
        [xgrid,ygrid,zgrid] = meshgrid(linspace(XL(1),XL(2),10),...
                    linspace(YL(1),YL(2),10),...
                    linspace(ZL(1),ZL(2),10));                     
        CO_b(:,:,1) = ones(10)*1; 
        CO_b(:,:,2) = ones(10)*1; 
        CO_b(:,:,3) = ones(10)*1; 
        surf(xgrid(:,:,1),ygrid(:,:,1),zgrid(:,:,1)*0+ZL_min, CO_b),
        shading interp

        
        % plot sperm tail=========================================================
        if i_nt>11
        for k=1:i_nt-11
            PlotTailCylinder_Movie_LabTrajProj(X(:,k),Y(:,k),Z(:,k),0.1,fig, XL(2),YL(2),ZL(1))
        end  
        end
        % Transparency parameter.
        FAlpha_mat_temp = [-0.9:0.1:-0.5 -0.45:0.05:-0.25 0]; %FAlpha_mat_temp = linspace(-1,0,11);
        FAlpha_mat = 10.^FAlpha_mat_temp;               
        i_alpha=0;
        for k=i_nt-10:i_nt
            i_alpha = i_alpha+1;
            PlotTailCylinder_Movie_LabTraj(X(:,k),Y(:,k),Z(:,k),FAlpha_mat(i_alpha),fig, XL(2),YL(2),ZL(1))
        end
        %
        % plot tail mid-point trajectory and projection=========================================================
        plot3(X(round(ns/2),1:i_nt),Y(round(ns/2),1:i_nt),Z(round(ns/2),1:i_nt),'color',[0.8500    0.3250    0.0980],'linewidth',1.5)    
        plot3(X(round(ns/2),1:i_nt),Y(round(ns/2),1:i_nt),Z(round(ns/2),1:i_nt)*0+ZL(1)+0.1,'color',[0.8500    0.3250    0.0980],'linewidth',0.1) 
        scatter3(X(round(ns/2),i_nt),Y(round(ns/2),i_nt),Z(round(ns/2),i_nt),60,[0.8500    0.3250    0.0980],'filled')    
        scatter3(X(round(ns/2),i_nt),Y(round(ns/2),i_nt),Z(round(ns/2),i_nt)*0+ZL(1)+0.1,30,[0.8500    0.3250    0.0980],'filled') 
            
        
        % Plot head==================================================================
        % Define head geometory.
        % neck (head centre) position to locate/ translate the sperm head
        B=[b1(:,i_nt),b2(:,i_nt),b3(:,i_nt)]; %orientation matrix (normalized) for sperm head
        neck=[X(2,i_nt),Y(2,i_nt),Z(2,i_nt)];          
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
        %
        source_ligth = [90 80];
        k = [0.9 1 1 5]; % [ka kd ks shine]. By default, k is [.55 .6 .4 10].
        s1= surfl(head_up_x,head_up_y,head_up_z,source_ligth,k );
        s2= surfl(head_down_x,head_down_y,head_down_z,source_ligth,k );
        s1.EdgeColor='none' ;
        s1.FaceColor=  [177 24 45]/255;  %[0.6 0 0];  
        s2.EdgeColor='none' ;
        s2.FaceColor=  [218 207 168]/255;  %[0 0.6 0];
        % head projection on XY plane.
        ss1= surf(head_up_x,head_up_y,head_up_z.*0+ZL(1) );
        ss2= surf(head_down_x,head_down_y,head_down_z.*0+ZL(1));
        cc = nan(size(s1.CData,1),size(s1.CData,2),3); cc(:,:,1)=0; cc(:,:,2)=0; cc(:,:,3)=0;ss1.CData=cc;
        cc = nan(size(s1.CData,1),size(s1.CData,2),3); cc(:,:,1)=0; cc(:,:,2)=0; cc(:,:,3)=0;ss2.CData=cc;             
        % head projection on XZ plane.
        ss1= surf(head_up_x,head_up_y.*0+YL(2),head_up_z);
        ss2= surf(head_down_x,head_down_y.*0+YL(2),head_down_z);
        cc = nan(size(s1.CData,1),size(s1.CData,2),3); cc(:,:,1)=0; cc(:,:,2)=0; cc(:,:,3)=0;ss1.CData=cc;
        cc = nan(size(s1.CData,1),size(s1.CData,2),3); cc(:,:,1)=0; cc(:,:,2)=0; cc(:,:,3)=0;ss2.CData=cc;             
        % head projection on YZ plane.
        ss1= surf(head_up_x.*0+XL(2),head_up_y,head_up_z);
        ss2= surf(head_down_x.*0+XL(2),head_down_y,head_down_z);
        cc = nan(size(s1.CData,1),size(s1.CData,2),3); cc(:,:,1)=0; cc(:,:,2)=0; cc(:,:,3)=0;ss1.CData=cc;
        cc = nan(size(s1.CData,1),size(s1.CData,2),3); cc(:,:,1)=0; cc(:,:,2)=0; cc(:,:,3)=0;ss2.CData=cc;             
           
        
        % Title: time.==================================================================      
        % (time in seconds)  
        dt = 1/90; % Each frame is taken every (1/90) seconds => 1.00 s corresponds to 90 video data frames (framerate was 90fps).
        inst = dt*(i_nt-1);
        format short
        inst = roundn(inst,-2);
        title(['Time: ',num2str(inst),'s'],'Fontname','Times New Roman','fontsize',fs,'fontweight','normal');
        
    
        xlim([XL(1),XL(2)]);ylim([YL(1),YL(2)]);zlim([ZL(1),ZL(2)]);
        xticks([-40:20:60]);yticks([-60:20:60]);zticks([-10:10:10])   
        ax = gca;ax.FontSize = fs; ax.TickLabelInterpreter = 'latex';
        xlh=xlabel('$x (\mu m)$ ','interpreter','latex','FontSize',1.5*fs);
        ylh=ylabel('$y (\mu m)$ ','interpreter','latex','FontSize',1.5*fs);
        zlh=zlabel('$z (\mu m)$ ','interpreter','latex','FontSize',1.5*fs,'Rotation',0);
        set(xlh,'position',[0 YL(1)-2 ZL(1)-5]);
        set(ylh,'position',[XL(1)-5 -10 ZL(1)-2]);
        set(zlh,'position',[XL(1)-8 YL(2)+5 -2]); % for no boundary %//////////////////////////////////////////////////////////////////
        %set(zlh,'position',[XL(1)-8 YL(2)+5 4]); % for 1 boundary %//////////////////////////////////////////////////////////////////
       
        
        l1=light;
        l1.Color = [1 1 1];
        l1.Style = 'local';  %'infinite';  %
        l1.Position = [-20 20 20];
        l2=light;
        l2.Color = [1 1 1];
        l2.Style = 'local';  %'infinite';  %
        l2.Position = [10 -10 20];
        l3=light;
        l3.Color = [1 1 1];
        l3.Style = 'local';  %'infinite';  %
        l3.Position = [40 -30 20];
        lighting gouraud   %phong  %flat  %
        material metal  %shiny     

        
        drawnow;
        frame(i_nt)=getframe(gcf);
        writeVideo(v,frame(i_nt));

   
    end
    
        
    end










%% Numerical reconstructed trace, with/ without a boundary.

if 1 %Video_num==1
    
    % Step 1: read numerical data (fine time discretization), and then get coarse-time interpolation and spatial smooth.
    load 01_lab_frame_raw_traces_2017_2018_rxmSSppxyz.mat; 
    nt = size(Z{sp},2);
    % (time in seconds)  
    dt=1/90; % Each frame is taken every (1/90) seconds => 1.00 s corresponds to 90 video data frames (framerate was 90fps).
    t_exp = 0:dt:(nt-1)*dt;
    for i_nt=1:nt
        arc_mat(:,i_nt) = arclength(X{sp}(:,i_nt), Y{sp}(:,i_nt), Z{sp}(:,i_nt));
    end
    arc_mean = mean(arc_mat(end,:)); % measured as 46.32.
    clearvars -except nt t_exp sp arc_mean 
    
    
    %load XNodes_sp8_NoBoundStokeslet.mat;  %//////////////////////////////////////////////////////////////////
    load XNodes_sp8_BoundBlakelet_H0.2.mat;    %//////////////////////////////////////////////////////////////////  
    X = Xs(Nhh+1:Nhh+Ns,:); Y = Xs(2*Nhh+Ns+1:2*Nhh+2*Ns,:); Z = Xs(3*Nhh+2*Ns+1:3*Nhh+3*Ns,:);
    nt_fine=size(X,2); %corresponding to the fine time discretization 'tRange' in the numerical reconstruction data, fine time for later flow field calcualtion.
    ns=size(X,1);
    head_tangent=head_tangent';
    head_normal=head_normal';
    head_binormal=nan(size(head_tangent)); % 3*nt
    for i_nt=1:nt_fine
        head_binormal(:,i_nt)=cross(head_tangent(:,i_nt),head_normal(:,i_nt));          
    end
    % Coarse-time interpolation.
    s_inp = linspace(0,1,ns);
    t_dim = t_nd/(2*pi)/Freq_WF;
    [Tcoarse,Scoarse] = meshgrid(t_exp,s_inp) ;% query points for coarser-time sperm tail coordinates.
    [Tfine,Sfine] = meshgrid(t_dim,s_inp) ; 
    X = interp2(Tfine,Sfine,X,Tcoarse,Scoarse,'spline'); % component 'x' of the tail coordinates, along arc length, at the instants of experimental recoarding 't_exp'.
    Y = interp2(Tfine,Sfine,Y,Tcoarse,Scoarse,'spline'); % component 'y' of the tail coordinates, along arc length, at the instants of experimental recoarding 't_exp'.
    Z = interp2(Tfine,Sfine,Z,Tcoarse,Scoarse,'spline'); % component 'z' of the tail coordinates, along arc length, at the instants of experimental recoarding 't_exp'.
    [Tcoarse,Scoarse] = meshgrid(t_exp,[1 2 3]) ;% query points for coarser-time sperm tail coordinates.
    [Tfine,Sfine] = meshgrid(t_dim,[1 2 3]) ; 
    head_tangent = interp2(Tfine,Sfine,head_tangent,Tcoarse,Scoarse,'spline');             
    head_normal  = interp2(Tfine,Sfine,head_normal,Tcoarse,Scoarse,'spline');             
    head_binormal = interp2(Tfine,Sfine,head_binormal,Tcoarse,Scoarse,'spline'); 
    b1 = head_tangent;
    b2 = head_normal;
    b3 = head_binormal;
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
    clearvars -except nt ns X Y Z b1 b2 b3 sp arc_mean ncycles 


    
    % Step 2: Generate sperm head: ellipsoid size semi-axes. 
    % The values can be modified, but their ratios should be unchanged.
    head_a1=1*arc_mean* 2.0/45; %
    head_a2=1*arc_mean* 1.6/45; %
    head_a3=1*arc_mean* 1.0/45; %
    %head:up+down£¬different colors,surf
    nth=20; 
    nphi=40; 
    a=1;
    M=nth*nphi;
    %head: up side
    [xup1,xup2,xup3]=GenerateSphereSurfaceForVisualisation_up(nth,nphi,a);
    [M1 M2]=size(xup1);
    xup1=xup1*head_a1; xup1=xup1-max(xup1(:));
    xup2=xup2*head_a2;
    xup3=xup3*head_a3;
    xup=[reshape(xup1,M,1);reshape(xup2,M,1);reshape(xup3,M,1)]; 
    %head: down side
    [xdown1,xdown2,xdown3]=GenerateSphereSurfaceForVisualisation_down(nth,nphi,a);
    xdown1=xdown1*head_a1; xdown1=xdown1-max(xdown1(:));
    xdown2=xdown2*head_a2;
    xdown3=xdown3*head_a3;
    xdown=[reshape(xdown1,M,1);reshape(xdown2,M,1);reshape(xdown3,M,1)]; 
    % locate the head axis for later location of head center in the lab frame.
    [~,np_head1] = max(xdown1(:));
    [~,np_head2] = min(xdown1(:));


    
        
    %% plot reconstructed lab trajectory, 'ns' space points and 'nt' time instants
    % Font options for printing figures
    fs=32;
    fn='times';
    fig = 3;
    
    figure(fig) 

    v=VideoWriter(['LabTraj_NumBound02'],'MPEG-4');  %//////////////////////////////////////////////////////////////////
    %v=VideoWriter(['LabTraj_NumNoBound'],'MPEG-4');  %//////////////////////////////////////////////////////////////////
    v.FrameRate=6;
    open(v);     

    % Locate the boundary, with the boundary distance to the initial head center position 0.2L/ infinite.                      
    % Define head geometory.
    % neck (head centre) position to locate/ translate the sperm head
    B=[b1(:,1),b2(:,1),b3(:,1)]; % (initial) orientation matrix (normalized) for sperm head
    neck=[X(1,1),Y(1,1),Z(1,1)];          
    % head surface rotation and translation
    % downside head
    xdownR=ApplyRotationMatrix(B,xdown);
    xdownT=TranslatePoints(xdownR,neck);
    [~,~,head_down_z]=ExtractComponents(xdownT);
    head_down_z_initial=reshape(head_down_z,M1,M2);
    % Locate the boundary.                      
    head_down_z_temp = head_down_z_initial(:);
    head_center_z = (head_down_z_temp(np_head1)+head_down_z_temp(np_head2))/2;
    % Set the boundary height.
    bound_height = 0.2+0.15;  %10;  %//////////////////////////////////////////////////////////////////
    % Note that for boundary height equal to 0.2L case, the visualized
    % boundary surface position is lower in Z coordinates, otherwise the
    % sperm body will contact with the boundary.
    ZL_min = head_center_z-bound_height; 
    
       
    X = X.*arc_mean; Y = Y.*arc_mean; Z = Z.*arc_mean;
    ZL_min = ZL_min*arc_mean;
    for i_nt= 11:nt-1  
        set(gcf,'color','w','units','normalized','outerposition',[0 0 1 1]) % make it white background
        clf;hold on;
        box on;axis equal;
        view(-35,20); % for sp#8          
        XL = [-0.9 1]*arc_mean;YL = [-0.8 0.5]*arc_mean;
        ZL = [ZL_min 0.5*arc_mean]; % for 1 boundary %//////////////////////////////////////////////////////////////////
        %ZL = [-0.3 0.25]*arc_mean; % for no boundary  %//////////////////////////////////////////////////////////////////
        
                       
        % Plot the boundary.==========================================================               
        [xgrid,ygrid,zgrid] = meshgrid(linspace(XL(1),XL(2),10),...
                    linspace(YL(1),YL(2),10),...
                    linspace(ZL(1),ZL(2),10));                     
        CO_b(:,:,1) = ones(10)*1; 
        CO_b(:,:,2) = ones(10)*1; 
        CO_b(:,:,3) = ones(10)*1; 
        surf(xgrid(:,:,1),ygrid(:,:,1),zgrid(:,:,1)*0+ZL_min, CO_b),
        shading interp

        
        % plot sperm tail=========================================================
        if i_nt>11
        for k=1:i_nt-11
            PlotTailCylinder_Movie_LabTrajProj(X(:,k),Y(:,k),Z(:,k),0.1,fig, XL(2),YL(2),ZL(1))
        end  
        end
        % Transparency parameter.
        FAlpha_mat_temp = [-0.9:0.1:-0.5 -0.45:0.05:-0.25 0]; %FAlpha_mat_temp = linspace(-1,0,11);
        FAlpha_mat = 10.^FAlpha_mat_temp;               
        i_alpha=0;
        for k=i_nt-10:i_nt
            i_alpha = i_alpha+1;
            PlotTailCylinder_Movie_LabTraj(X(:,k),Y(:,k),Z(:,k),FAlpha_mat(i_alpha),fig, XL(2),YL(2),ZL(1))
        end
        %
        % plot tail mid-point trajectory and projection=========================================================
        plot3(X(round(ns/2),1:i_nt),Y(round(ns/2),1:i_nt),Z(round(ns/2),1:i_nt),'color',[0.8500    0.3250    0.0980],'linewidth',1.5)    
        plot3(X(round(ns/2),1:i_nt),Y(round(ns/2),1:i_nt),Z(round(ns/2),1:i_nt)*0+ZL(1)+0.1,'color',[0.8500    0.3250    0.0980],'linewidth',0.1) 
        scatter3(X(round(ns/2),i_nt),Y(round(ns/2),i_nt),Z(round(ns/2),i_nt),60,[0.8500    0.3250    0.0980],'filled')    
        scatter3(X(round(ns/2),i_nt),Y(round(ns/2),i_nt),Z(round(ns/2),i_nt)*0+ZL(1)+0.1,30,[0.8500    0.3250    0.0980],'filled') 
            
        
        % Plot head==================================================================
        % Define head geometory.
        % neck (head centre) position to locate/ translate the sperm head
        B=[b1(:,i_nt),b2(:,i_nt),b3(:,i_nt)]; %orientation matrix (normalized) for sperm head
        neck=[X(2,i_nt),Y(2,i_nt),Z(2,i_nt)];          
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
        %
        source_ligth = [90 80];
        k = [0.9 1 1 5]; % [ka kd ks shine]. By default, k is [.55 .6 .4 10].
        s1= surfl(head_up_x,head_up_y,head_up_z,source_ligth,k );
        s2= surfl(head_down_x,head_down_y,head_down_z,source_ligth,k );
        s1.EdgeColor='none' ;
        s1.FaceColor=  [177 24 45]/255;  %[0.6 0 0];  
        s2.EdgeColor='none' ;
        s2.FaceColor=  [218 207 168]/255;  %[0 0.6 0];
        % head projection on XY plane.
        ss1= surf(head_up_x,head_up_y,head_up_z.*0+ZL(1) );
        ss2= surf(head_down_x,head_down_y,head_down_z.*0+ZL(1));
        cc = nan(size(s1.CData,1),size(s1.CData,2),3); cc(:,:,1)=0; cc(:,:,2)=0; cc(:,:,3)=0;ss1.CData=cc;
        cc = nan(size(s1.CData,1),size(s1.CData,2),3); cc(:,:,1)=0; cc(:,:,2)=0; cc(:,:,3)=0;ss2.CData=cc;             
        % head projection on XZ plane.
        ss1= surf(head_up_x,head_up_y.*0+YL(2),head_up_z);
        ss2= surf(head_down_x,head_down_y.*0+YL(2),head_down_z);
        cc = nan(size(s1.CData,1),size(s1.CData,2),3); cc(:,:,1)=0; cc(:,:,2)=0; cc(:,:,3)=0;ss1.CData=cc;
        cc = nan(size(s1.CData,1),size(s1.CData,2),3); cc(:,:,1)=0; cc(:,:,2)=0; cc(:,:,3)=0;ss2.CData=cc;             
        % head projection on YZ plane.
        ss1= surf(head_up_x.*0+XL(2),head_up_y,head_up_z);
        ss2= surf(head_down_x.*0+XL(2),head_down_y,head_down_z);
        cc = nan(size(s1.CData,1),size(s1.CData,2),3); cc(:,:,1)=0; cc(:,:,2)=0; cc(:,:,3)=0;ss1.CData=cc;
        cc = nan(size(s1.CData,1),size(s1.CData,2),3); cc(:,:,1)=0; cc(:,:,2)=0; cc(:,:,3)=0;ss2.CData=cc;             
           
        
        % Title: time.==================================================================      
        % (time in seconds)  
        dt = 1/90; % Each frame is taken every (1/90) seconds => 1.00 s corresponds to 90 video data frames (framerate was 90fps).
        inst = dt*(i_nt-1);
        format short
        inst = roundn(inst,-2);
        title(['Time: ',num2str(inst),'s'],'Fontname','Times New Roman','fontsize',fs,'fontweight','normal');
        
    
        xlim([XL(1),XL(2)]);ylim([YL(1),YL(2)]);zlim([ZL(1),ZL(2)]);
        xticks([-40:20:60]);yticks([-60:20:60]);zticks([-10:10:10])   
        ax = gca;ax.FontSize = fs; ax.TickLabelInterpreter = 'latex';
        xlh=xlabel('$x (\mu m)$ ','interpreter','latex','FontSize',1.5*fs);
        ylh=ylabel('$y (\mu m)$ ','interpreter','latex','FontSize',1.5*fs);
        zlh=zlabel('$z (\mu m)$ ','interpreter','latex','FontSize',1.5*fs,'Rotation',0);
        set(xlh,'position',[0 YL(1)-2 ZL(1)-5]);
        set(ylh,'position',[XL(1)-5 -10 ZL(1)-2]);
        %set(zlh,'position',[XL(1)-8 YL(2)+5 -2]); % for no boundary %//////////////////////////////////////////////////////////////////
        set(zlh,'position',[XL(1)-8 YL(2)+5 4]); % for 1 boundary %//////////////////////////////////////////////////////////////////
       
        
        l1=light;
        l1.Color = [1 1 1];
        l1.Style = 'local';  %'infinite';  %
        l1.Position = [-20 20 20];
        l2=light;
        l2.Color = [1 1 1];
        l2.Style = 'local';  %'infinite';  %
        l2.Position = [10 -10 20];
        l3=light;
        l3.Color = [1 1 1];
        l3.Style = 'local';  %'infinite';  %
        l3.Position = [40 -30 20];
        lighting gouraud   %phong  %flat  %
        material metal  %shiny     

        
        drawnow;
        frame(i_nt)=getframe(gcf);
        writeVideo(v,frame(i_nt));

   
    end
    
        
    end

end
