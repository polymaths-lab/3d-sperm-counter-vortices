function Generate_FlowStreamlines(Xg,Yg,Zg,Ug_mean,Vg_mean,Wg_mean,sx_3D,sy_3D,sz_3D, vel)


% Locate each vertices along the 3D streamlines.        
vert = stream3(Xg,Yg,Zg,Ug_mean,Vg_mean,Wg_mean,sx_3D,sy_3D,sz_3D);


% Get the velocity values along streamlines.
% 'Vetices' is a cell array and each cell in this case has size '?*3', which gives the [X Y Z] coordinates along each streamline.
% That's why we need to interpolate the velocity field to get the velocity
% at these vertices.
vert_mat = [0 0 0];
nsl = length(vert);  %number of streamlines
nvert_sl = nan(nsl,1); % record the number of vertices for each streamline.
for i_sl = 1:nsl
    vert_temp = vert{i_sl};
    nvert_sl(i_sl) = size(vert_temp,1);
    vert_mat = [vert_mat; vert_temp];
end
vert_mat = vert_mat(2:end,:);


% Map streamline colors according to velocity values.
Vq = interp3(Xg,Yg,Zg,vel,vert_mat(:,1),vert_mat(:,2),vert_mat(:,3));
nv = 0;
for i_sl = 1:nsl
    nv0 = nv+1;
    nv = nv+nvert_sl(i_sl);
    vert_temp = { vert{i_sl} };
    %h = streamtube(vert_temp, 0.01); %near full flow for idealized model--> plot.
    h = streamtube(vert_temp, 0.03); %near full flow for idealized model--> video.
    %h = streamtube(vert_temp, 0.05); %near cross-sectional flow for idealized model.
    %h = streamtube(vert_temp, 0.3); %far flow for idealized model.
    h.CData = repmat(Vq(nv0:nv),1,21);   
end

colormap(gca,'parula')
    