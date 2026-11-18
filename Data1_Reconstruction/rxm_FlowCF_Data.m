function rxm_FlowCF_Data

% This function calculates the flow field around a sperm, relative to the
% comovinf frame of reference, where the flow field translates but not 
% rotates with the sperm head, and the X-axis is aligned with the  
% time-averaged flagellar position

% By adjusting the parameters, the calculated flow data can be for near or
% far field, with or without a boundary, coarse or fine time points.

% Note: the sperm head center is at the origin of the comoving-frame field.

clc

if 1 % Get the flow data.
% near-boundary sperm.
sp=8;

load FreeSperm_Frequency.mat;
Freq_WF = HF_Freq(sp);
Freq_rev = CF_Freq(sp);


%% Step1: define the flow field size and discretization (near/ far field).

FieldRange = 'n'; % near-field
%FieldRange = 'f'; % far-field

if FieldRange == 'n'
    xl_field=[-1.5 0.5]; yl_field=[-1 1]; zl_field=[-1 1];
    Nx_field=21; Ny_field=21; Nz_field=21;
elseif FieldRange == 'f'
    % For 'No Bound' case
    %xl_field=[-11 10]; yl_field=[-10 10]; zl_field=[-10 10];
    %Nx_field=43; Ny_field=41; Nz_field=41;
    % For '1 Bound' case
    xl_field=[-11 10]; yl_field=[-10 10]; zl_field=[-1 10];
    Nx_field=43; Ny_field=41; Nz_field=23;
end

% 3D grid space of the local flow field.
xg=linspace(xl_field(1),xl_field(2),Nx_field);
yg=linspace(yl_field(1),yl_field(2),Ny_field); 
zg=linspace(zl_field(1),zl_field(2),Nz_field); 
[Xg_raw,Yg_raw,Zg_raw]=meshgrid(xg,yg,zg);



%% Step2: set the position of the flow field relative to the sperm body--comoving frame (coarse/ fine time).

% First, align the field box to the orientation, where its X axis is
% parallel to the average flagellar direction in 'comoving' frame.
load XNodes_sp8_NoBoundStokeslet_Tfine.mat; %
%load XNodes_sp8_BoundBlakelet_H0.5.mat; %////////////////////////////////////////////////////////////////////////////////////////////////////////
Q=Nhh+Ns;
Xtail_0 = Xs(Nhh+1:Q,:); 
Ytail_0 = Xs(Q+Nhh+1:2*Q,:); 
Ztail_0 = Xs(2*Q+Nhh+1:3*Q,:); 
%
% 以下对鞭毛形状光滑处理的代码是新添的，想看一下该处理能否获得较规律的流场PCA结果。在此之前的计算步骤没有这一步是因为，发现未经过光滑处理的波形重构所得的精子轨迹，相比于光滑处理的波形，呈现出更好的重构精度。
q=nhh+ns;
xtail_0 = xs(nhh+1:q,:); 
ytail_0 = xs(q+nhh+1:2*q,:); 
ztail_0 = xs(2*q+nhh+1:3*q,:); 
% Spatial smooth of the flagellar shape, to remove the noise arising from the raw experimental data being reconstructed.
    ppx = 0.999; % 'ppx=1' means interpolation without cubic spline.
    ppy = 0.999;
    ppz = 0.999;  
    nt_XNodes = size(Xtail_0,2);
    for i_nt = 1:nt_XNodes
        S_temp = arclength(Xtail_0(:,i_nt), Ytail_0(:,i_nt), Ztail_0(:,i_nt));
        Xtail_0(:,i_nt) = fnval(csaps(S_temp,Xtail_0(:,i_nt),ppx),S_temp);
        Ytail_0(:,i_nt) = fnval(csaps(S_temp,Ytail_0(:,i_nt),ppy),S_temp);
        Ztail_0(:,i_nt) = fnval(csaps(S_temp,Ztail_0(:,i_nt),ppz),S_temp);
        s_temp = arclength(xtail_0(:,i_nt), ytail_0(:,i_nt), ztail_0(:,i_nt));
        xtail_0(:,i_nt) = fnval(csaps(s_temp,xtail_0(:,i_nt),ppx),s_temp);
        ytail_0(:,i_nt) = fnval(csaps(s_temp,ytail_0(:,i_nt),ppy),s_temp);
        ztail_0(:,i_nt) = fnval(csaps(s_temp,ztail_0(:,i_nt),ppz),s_temp);
    end
%
Xtail_1=Xtail_0-repmat(Xtail_0(1,:),Ns,1);    
Ytail_1=Ytail_0-repmat(Ytail_0(1,:),Ns,1);    
Ztail_1=Ztail_0-repmat(Ztail_0(1,:),Ns,1);    
xtail_mean= mean(Xtail_1,2); %Ns*1
ytail_mean= mean(Ytail_1,2);
ztail_mean= mean(Ztail_1,2);

dir_raw=[1 0 0];
dir_new=[xtail_mean(1)-xtail_mean(end),ytail_mean(1)-ytail_mean(end),ztail_mean(1)-ztail_mean(end)];

[Xg_align_temp,Yg_align_temp,Zg_align_temp] = rxm_dir_align(dir_raw,dir_new,Xg_raw(:),Yg_raw(:),Zg_raw(:));
Xg_align = reshape(Xg_align_temp,size(Xg_raw,1),size(Xg_raw,2),size(Xg_raw,3));
Yg_align = reshape(Yg_align_temp,size(Xg_raw,1),size(Xg_raw,2),size(Xg_raw,3));
Zg_align = reshape(Zg_align_temp,size(Xg_raw,1),size(Xg_raw,2),size(Xg_raw,3));




% Second, translate the aligned field box to make it moves with sperm head
% center.
nBC_raw = t_nd(end)/(2*pi); %t_nd: non-dimensioanl time.
nt_raw = length(t_nd);
t_dim = t_nd/(2*pi)/Freq_WF; %dimensional time, in units of 's'.
np_beat_raw = nt_raw/nBC_raw;
% Two purposes here:
%    Select an appropriate time range, i.e. complete trace revolution periods.
%    Time interpolation, because in some cases, the flow data has coarser time points than the cell mobility data.
if 0 % Coarse-time flow data for plots, such as average flow field, which does not require fine temporal discretization, so using coarse time discretization.
    np_beat_temp=20;
    dnp = round(np_beat_raw/np_beat_temp);
    % Find the time range of complete comoving trace revolutions.
    nrev_temp = floor(Freq_rev*t_dim(end));
else % Fine-time flow data for PCA and video, which require fine temporal discretization.
    np_beat_temp=np_beat_raw;
    dnp = round(np_beat_raw/np_beat_temp);
    if 1 % Find the time range of complete comoving trace revolutions, for PCA.
        nrev_temp = floor(Freq_rev*t_dim(end));
    else % Find the time range for just 1 comoving trace revolution, for video.
        nrev_temp = 1;
    end
end
[~,nt_temp] = min(abs(Freq_rev*t_dim-nrev_temp));
t_nd_flow = t_nd(1:dnp:nt_temp);
nt_flow = length(t_nd_flow); 
xfield = [Xg_align(:); Yg_align(:); Zg_align(:)];
Xfield = nan(length(xfield),nt_flow);
for i_nt=1:nt_flow  
    x0=head_x0(1+(i_nt-1)*dnp,:);
    xfield_temp=TranslatePoints(xfield,x0);
    Xfield(:,i_nt)=xfield_temp;
end




%% Step3: Calculate instantaneous flow field (with/ without a boundary).

% numerical parameters
    epsilon=0.25/45;
    domain='i';  %Stokeslet %////////////////////////////////////////////////////////////////////////////////////////////////////////
    %domain='h';    %Blakelet %////////////////////////////////////////////////////////////////////////////////////////////////////////
    blockSize=0.2;

    Ug = nan(length(Xg_align(:)),nt_flow); 
    Vg = nan(length(Yg_align(:)),nt_flow); 
    Wg = nan(length(Zg_align(:)),nt_flow);
    for i_nt=1:nt_flow 
        i_nt
        % Assemble allocation points: sperm and boundary
        %x=xs(:,1+(i_nt-1)*dnp);%   x=MergeVectorGrids(xs(:,1+(i_nt-1)*dnp),xb); <-- 'xb=[]'.
        %X=Xs(:,1+(i_nt-1)*dnp);%   X=MergeVectorGrids(Xs(:,1+(i_nt-1)*dnp),Xb); <-- 'xb=[]'.
        %
        % 此处代码是新改的。目的在于计算光滑处理后的鞭毛所产生的流场，也许这样的流场能够得到更规律的PCA结果。在此之前的计算步骤没有这一步，是因为精子轨迹运动的重构精度显示，未光滑处理的波形相比于光滑处理的波形显示出更高的动力重构精度。
        x=[xs(1:nhh,1+(i_nt-1)*dnp); xtail_0(:,1+(i_nt-1)*dnp); xs(q+1:q+nhh,1+(i_nt-1)*dnp);  ytail_0(:,1+(i_nt-1)*dnp);  xs(2*q+1:2*q+nhh,1+(i_nt-1)*dnp); ztail_0(:,1+(i_nt-1)*dnp)];
        X=[Xs(1:Nhh,1+(i_nt-1)*dnp); Xtail_0(:,1+(i_nt-1)*dnp);  Xs(Q+1:Q+Nhh,1+(i_nt-1)*dnp); Ytail_0(:,1+(i_nt-1)*dnp);  Xs(2*Q+1:2*Q+Nhh,1+(i_nt-1)*dnp); Ztail_0(:,1+(i_nt-1)*dnp)];        
        %
        % evaluate velocity field on grid and then save the data
        u=EvaluateVelocityFromForce(Xfield(:,i_nt),X,x,f(:,1+(i_nt-1)*dnp),epsilon,domain,blockSize);
        [ug,vg,wg]=ExtractComponents(u);

        Ug(:,i_nt) = ug; Vg(:,i_nt) = vg; Wg(:,i_nt) = wg;
    end    
   
    
 
%////////////////////////////////////////////////////////////////////////////////////////////////////////
%save('FarFlow_CF_sp8_NoBound_FineTimeForVideo.mat',...
%save('NearFlow_CF_sp8_BoundH0.2.mat',...
%save('NearFlow_CF_sp8_NoBound.mat',...
%save('FarFlow_CF_sp8_NoBound.mat',...
save('NearFlow_CF_sp8_NoBound_Tfine.mat',...
    'xl_field','yl_field','zl_field','Nx_field','Ny_field','Nz_field',...
 'Ug','Vg','Wg','t_nd_flow','Freq_WF','Freq_rev');


end


%% Get flow PCA.

if 1
    load NearFlow_CF_sp8_NoBound_Tfine.mat;

    %The computational cost exceeds the allowed memory, so coarse
    %interpolation is needed here.
    Ug = Ug(:,1:2:end);
    Vg = Vg(:,1:2:end);
    Wg = Wg(:,1:2:end);
    
    
    np = size(Ug,1);
    [PCA_mode,Time_coef, nmode,cumulative_variance]=Extract_flow_PCA(Ug,Vg,Wg); 

    PCA_mode_Ug = PCA_mode(1:np,:); 
    PCA_mode_Vg = PCA_mode(np+1:2*np,:); 
    PCA_mode_Wg = PCA_mode(2*np+1:3*np,:);
    
    nmode
    
save(['NearFlow_CF_sp8_NoBound_Tfine_PCA_Accuracy99.99.mat'],'PCA_mode_Ug','PCA_mode_Vg','PCA_mode_Wg','Time_coef','nmode','cumulative_variance');  

end


