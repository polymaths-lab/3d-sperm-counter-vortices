function [xs,Xs] = GetQuadForceGrid(B,x0,Xh,xh,Xt,xt,tspan,t)


    nt = size(Xt,2); % time interval points
    Ns = size(Xt,1)/3; % 'Ns' for fine grid - position only
    ns = size(xt,1)/3; % 'ns' for coarse grid - position and velocity  


    dt = 0.001; % 'dt' for calculating velocity
    t_inp = tspan; 
    S_inp = linspace(0,1,Ns);
    s_inp = linspace(0,1,ns);
    [T,S] = meshgrid(t_inp,S_inp) ;
    [Tq,Sq] = meshgrid(t,S_inp) ; % query points for 'Xt'
    [Tt,s] = meshgrid(t_inp,s_inp) ;
    %[Ttq,sq] = meshgrid([t-dt/2 t t+dt/2],s_inp) ; % query points for 'xt'
    [Ttq,sq] = meshgrid(t,s_inp) ; % query points for 'xt'
    
    

    % Interpolation:
    Xt1 = Xt(1:Ns,:); Xt2 = Xt(1+Ns:2*Ns,:); Xt3 = Xt(1+2*Ns:3*Ns,:); %sperm tail coordinates (fine grid)
    Xt1q = interp2(T,S,Xt1,Tq,Sq); % component 'x' of the tail coordinates, along arc length(S,Sq) at the instant 'Tq' 
    Xt2q = interp2(T,S,Xt2,Tq,Sq); % component 'y' of the tail coordinates, along arc length(S,Sq) at the instant 'Tq'
    Xt3q = interp2(T,S,Xt3,Tq,Sq); % component 'z' of the tail coordinates, along arc length(S,Sq) at the instant 'Tq'
    xt1 = xt(1:ns,:); xt2 = xt(1+ns:2*ns,:); xt3 = xt(1+2*ns:3*ns,:);%sperm tail coordinates (coarse grid)
    xt1q = interp2(Tt,s,xt1,Ttq,sq);%,'spline'); 
    xt2q = interp2(Tt,s,xt2,Ttq,sq);%,'spline'); 
    xt3q = interp2(Tt,s,xt3,Ttq,sq);%,'spline'); 
% IMPORTANT!
% 'spline' method is important for 'interp2' when it comes with extrapval
% data; otherwise, NAN value will be obtained 
    %xt1qq = interp2(Tt,s,xt1,Ttq,sq,'spline'); 
    %xt2qq = interp2(Tt,s,xt2,Ttq,sq,'spline'); 
    %xt3qq = interp2(Tt,s,xt3,Ttq,sq,'spline'); 
    %xt1q = xt1qq(:,2); %sperm tail coordinates (coarse grid)
    %xt2q = xt2qq(:,2);
    %xt3q = xt3qq(:,2);
    

% coordinates of the sperm head, including [xh,vh,Xh]======================================================== 
    [Xh1,Xh2,Xh3]=ExtractComponents(Xh);
    [xh1,xh2,xh3]=ExtractComponents(xh);

% Combine the sperm head and sperm tail information===============================================================    
% Since the position of sperm head in the body frame is always fixed, we only need to query the changing positions of tail, such as 'xt1q','vt1q'.  
    xi = [xh1;xt1q;xh2;xt2q;xh3;xt3q];
    Xi = [Xh1;Xt1q;Xh2;Xt2q;Xh3;Xt3q];


    xx0=ApplyRotationMatrix(B,xi); % x-x0
    xs=TranslatePoints(xx0,x0);
    Xx0=ApplyRotationMatrix(B,Xi); % x-x0
    Xs=TranslatePoints(Xx0,x0);