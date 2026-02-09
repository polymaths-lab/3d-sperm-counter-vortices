function rxm_DevDis_VarySperm_Plot

clc

%% Read the data of the sperm displacement relative deviation between with-boundary and no-boundary cases.

load XNodes_sp8_NoBoundStokeslet.mat; 
dis_temp = arclength(head_x0(:,1),head_x0(:,2),head_x0(:,3)); 
dis_NB = dis_temp(end);
clearvars -except dis_NB
load XNodes_sp8_BoundBlakelet_H0.2.mat;
dis_temp = arclength(head_x0(:,1),head_x0(:,2),head_x0(:,3)); 
dis_BB = dis_temp(end);
Dis_RelDev_exp_raw = (dis_BB-dis_NB)/dis_NB
clearvars -except Dis_RelDev_exp_raw


load DisDev_VarySperm_31BC.mat;
Dis_RelDev_virtual_raw = Dis_RelDev


%% Barplot for the displacement relative deviation.

scal = 30;
Dis_RelDev_exp = Dis_RelDev_exp_raw*scal; 
Dis_RelDev_virtual = Dis_RelDev_virtual_raw*scal; 
clrs_min = min([Dis_RelDev_exp; Dis_RelDev_virtual(:)]);
clrs_max = max([Dis_RelDev_exp; Dis_RelDev_virtual(:)]);
nclrs = 100; 
clrs_map = parula(nclrs);


fs = 40;
figure(1)
clf; hold on; axis equal; set(gcf,'color','w')
view([-50 30]);  %view(3); 
box on; grid on;
ax = gca; ax.FontSize = fs;  ax.TickLabelInterpreter = 'latex';



%% Plot for all the virtual sperm cases.

n_col = size(Dis_RelDev_virtual,2); %n_alpha
n_row = size(Dis_RelDev_virtual,1); %n_head
bar_virtual = bar3(Dis_RelDev_virtual,0.6);  %1*n_col. All elemets are surfaces.

for i_col = 1:n_col
    
    bar_surf_virtual = bar_virtual(i_col); % (n_row*6) *4: 6 surfaces for each bar, and 4 vertex for each surface.
    bar_surf_x = bar_surf_virtual.XData; 
    bar_surf_y = bar_surf_virtual.YData; 
    bar_surf_z = bar_surf_virtual.ZData; 
    bar_surf_color = nan( (n_row*6), 4, 3);
    
    for i_row = 1:n_row  
        
        bar_surf_z_temp0 = bar_surf_z((i_row-1)*6+1:i_row*6,:); 
        bar_surf_z_temp = bar_surf_z_temp0;
        ind_NAN = ~isnan( bar_surf_z_temp ); % for NaN values, ind_NAN=0; otherwise =1.
        bar_surf_z_temp(find(~ind_NAN)) = 0;
        bar_surf_z_temp0( find(bar_surf_z_temp) ) = abs( Dis_RelDev_virtual(i_row,i_col) );
        bar_surf_z((i_row-1)*6+1:i_row*6,:) = bar_surf_z_temp0; 
                
        clrs_ind = ( Dis_RelDev_virtual(i_row,i_col)-clrs_min )/ (clrs_max-clrs_min);
        clrs_ind = round(clrs_ind*nclrs);
        if clrs_ind==0
            clrs_ind=1;
        end
        bar_surf_color((i_row-1)*6+1:i_row*6,:,1) = repmat(clrs_map(clrs_ind,1),6,4);
        bar_surf_color((i_row-1)*6+1:i_row*6,:,2) = repmat(clrs_map(clrs_ind,2),6,4);
        bar_surf_color((i_row-1)*6+1:i_row*6,:,3) = repmat(clrs_map(clrs_ind,3),6,4);
       
        
        
        if i_row==1 && i_col==1 
            % Extract the minimum X and Y limits.
            xrange_temp = bar_surf_x((i_row-1)*6+1:i_row*6,:); 
            yrange_temp = bar_surf_y((i_row-1)*6+1:i_row*6,:); 
            XL_min = min(xrange_temp(:));
            YL_min = min(yrange_temp(:));
        elseif i_row==1 && i_col==4
            % For the location of the bar of our observation case.
            bar_surf_x_temp1 = bar_surf_x((i_row-1)*6+1:i_row*6,:);
            bar_surf_y_temp1 = bar_surf_y((i_row-1)*6+1:i_row*6,:);
        elseif i_row==1 && i_col==5
            % For the location of the bar of our observation case.
            bar_surf_x_temp2 = bar_surf_x((i_row-1)*6+1:i_row*6,:);
        elseif i_row==n_row && i_col==n_col
            % Extract the maximum X and Y limits.
            xrange_temp = bar_surf_x((i_row-1)*6+1:i_row*6,:); 
            yrange_temp = bar_surf_y((i_row-1)*6+1:i_row*6,:); 
            XL_max = max(xrange_temp(:));
            YL_max = max(yrange_temp(:));
        end
    end
    
    bar_surf_virtual.ZData = bar_surf_z;
    bar_surf_virtual.CData = bar_surf_color;
    bar_surf_virtual.EdgeColor = [1 1 1]*0.8;%bar_surf_color(1,1,:);
    bar_surf_virtual.FaceAlpha = 0.3;
end


%% Plot the observed sperm case.

bar_surf_exp = bar3(Dis_RelDev_exp,0.1);  %1*1. All elemets are surfaces.

% Adjust the bar position: [X Y].
bar_surf_exp.XData = bar_surf_x_temp1 + (bar_surf_x_temp2-bar_surf_x_temp1).*0.3;
bar_surf_exp.YData = bar_surf_y_temp1;
bar_x_min = min(bar_surf_exp.XData(:));  bar_x_max = max(bar_surf_exp.XData(:)); bar_x_mean = (bar_x_max+bar_x_min)/2;
bar_y_min = min(bar_surf_exp.YData(:));  bar_y_max = max(bar_surf_exp.YData(:)); bar_y_mean = (bar_y_max+bar_y_min)/2;
bar_width = 0.3;
bar_surf_exp.XData(find(bar_surf_exp.XData==bar_x_min)) = bar_x_mean-bar_width/2;
bar_surf_exp.XData(find(bar_surf_exp.XData==bar_x_max)) = bar_x_mean+bar_width/2;
bar_surf_exp.YData(find(bar_surf_exp.YData==bar_y_min)) = bar_y_mean-bar_width/2;
bar_surf_exp.YData(find(bar_surf_exp.YData==bar_y_max)) = bar_y_mean+bar_width/2;

% Adjust the bar height: Z.
bar_surf_z_temp0 = bar_surf_exp.ZData; 
bar_surf_z_temp = bar_surf_z_temp0;
ind_NAN = ~isnan( bar_surf_z_temp ); % for NaN values, ind_NAN=0; otherwise =1.
bar_surf_z_temp(find(~ind_NAN)) = 0;
bar_surf_z_temp0( find(bar_surf_z_temp) ) = abs( Dis_RelDev_exp );
bar_surf_z = bar_surf_z_temp0; 
bar_surf_exp.ZData = bar_surf_z;

% Adjust the bar color.
bar_surf_color = nan(6,4,3);  % each bar has 6 faces; each face has 4 vertex; each color is defined by 3 numbers.
clrs_ind = ( Dis_RelDev_exp-clrs_min )/ (clrs_max-clrs_min);
clrs_ind = round(clrs_ind*nclrs);
if clrs_ind==0  
    clrs_ind=1;
end
bar_surf_color(:,:,1) = repmat(clrs_map(clrs_ind,1),6,4);
bar_surf_color(:,:,2) = repmat(clrs_map(clrs_ind,2),6,4);
bar_surf_color(:,:,3) = repmat(clrs_map(clrs_ind,3),6,4);
bar_surf_exp.CData = bar_surf_color;                        

% Other bar setting.
bar_surf_exp.EdgeColor = [1 1 1]*0;
bar_surf_exp.LineWidth=1;
bar_surf_exp.FaceAlpha = 1;


%% Coordinate system outlook setting.

XL=[XL_min XL_max]; YL=[YL_min YL_max]; ZL=[0 13]; 
xlim([XL(1),XL(2)]);ylim([YL(1),YL(2)]); zlim([ZL(1),ZL(2)]); 
xticks([1:2:5]); yticks([1:4]); zticks([0:0.2:0.4].*scal); 
xticklabels({'0','0.4','0.8'}); yticklabels({'1','5','10','15'}); zticklabels({'0','20','40'}); 
xlh=xlabel('$\alpha$','interpreter','latex','FontSize',fs+10);
ylh=ylabel('Head scaling factor','interpreter','latex','FontSize',fs);
zlh=zlabel('$|\Delta D| \ (100\%)$','interpreter','latex','FontSize',fs,'Rotation',90);
%zlh=zlabel('$\Delta D (100\%)$','interpreter','latex','FontSize',fs,'Rotation',90);
 

colormap(parula)
c= colorbar; caxis([clrs_min clrs_max]./scal*100); 
set(c,'position',[0.3,0.15,0.01,0.7]); %c.Location = 'west'; 
c.Label.String = '\Delta D'; c.TickLabelInterpreter='latex'; 
c.Ticks=[-40 -20 0]; c.TickLabels={'-40','-20','0'}; 
c.FontName='times'; c.FontSize=fs;



