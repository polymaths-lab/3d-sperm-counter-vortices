function dz=SolveSwimmingProblemWithBoundary_WF(z,Xh,xh,vh,Xt,xt,boundary,t,tspan,epsilon,domain,blockSize,varargin)

    % solves force-free swimming problem for translational and angular velocity U, Om
    % in the presence of a stationary rigid boundary (e.g. finite plane wall)
    %
    % input: z       - position/orientation of swimmer
    %                  dz(1:3) = x0      - origin of swimmer
    %                  dz(4:6) = b1      - first basis vector of swimmer frame
    %                  dz(7:9) = b2      - second basis vector of swimmer frame
    %        swimmer - structure describing how to construct swimmer
    %        xb      - discretisation of stationary boundary - force points
    %        Xb      - discretisation of stationary boundary - quadrature points
    %        t       - time (scalar)
    %
    % variables: U      - translational velocity
    %            Om     - angular velocity
    %
    % output: dz        - rate of change of position-orientation (x,b1,b2)
    %                     first three compts are dx0/dt, 
    %                     next three are db1/dt,
    %                     next three are db2/dt
    %                   - if varargin{1}='f' then includes force components

    x0=z(1:3);
    b1=z(4:6);
    b2=z(7:9);
    b3=cross(b1,b2);
    B=[b1(:) b2(:) b3(:)];


%% Position of the sperm head and tail in the body frame, including [xi,vi,Xi]
% interpolation of the tail waveform coordinates, including [xt,vt,Xt]=========================================================
    %nt = size(Xt,2); % time interval points
    Ns = size(Xt,1)/3; % 'Ns' for fine grid - position only
    ns = size(xt,1)/3; % 'ns' for coarse grid - position and velocity  

    dt = 0.001; % 'dt' for calculating velocity
    t/(2*pi)
    t_inp = tspan; 
    S_inp = linspace(0,1,Ns);
    s_inp = linspace(0,1,ns);
    [T,S] = meshgrid(t_inp,S_inp) ;
    [Tq,Sq] = meshgrid(t,S_inp) ; % query points for 'Xt'
    [Tt,s] = meshgrid(t_inp,s_inp) ;
    [Ttq,sq] = meshgrid([t-dt/2 t t+dt/2],s_inp) ; % query points for 'xt'
    
    % Interpolation:
    Xt1 = Xt(1:Ns,:); Xt2 = Xt(1+Ns:2*Ns,:); Xt3 = Xt(1+2*Ns:3*Ns,:); %sperm tail coordinates (fine grid)
    Xt1q = interp2(T,S,Xt1,Tq,Sq); % component 'x' of the tail coordinates, along arc length(S,Sq) at the instant 'Tq' 
    Xt2q = interp2(T,S,Xt2,Tq,Sq); % component 'y' of the tail coordinates, along arc length(S,Sq) at the instant 'Tq'
    Xt3q = interp2(T,S,Xt3,Tq,Sq); % component 'z' of the tail coordinates, along arc length(S,Sq) at the instant 'Tq'
    xt1 = xt(1:ns,:); xt2 = xt(1+ns:2*ns,:); xt3 = xt(1+2*ns:3*ns,:);%sperm tail coordinates (coarse grid)
    % 'spline' method is important for 'interp2' when it comes with extrapval data; otherwise, NAN value will be obtained. 
    xt1qq = interp2(Tt,s,xt1,Ttq,sq,'spline'); 
    xt2qq = interp2(Tt,s,xt2,Ttq,sq,'spline'); 
    xt3qq = interp2(Tt,s,xt3,Ttq,sq,'spline'); 
    xt1q = xt1qq(:,2); %sperm tail coordinates (coarse grid)
    xt2q = xt2qq(:,2);
    xt3q = xt3qq(:,2);
    vt1q = (xt1qq(:,3)-xt1qq(:,1))/dt; %sperm tail velocity (coarse grid)
    vt2q = (xt2qq(:,3)-xt2qq(:,1))/dt;
    vt3q = (xt3qq(:,3)-xt3qq(:,1))/dt;

    % coordinates of the sperm head, including [xh,vh,Xh]========================================================
    [Xh1,Xh2,Xh3]=ExtractComponents(Xh);
    [xh1,xh2,xh3]=ExtractComponents(xh);
    [vh1,vh2,vh3]=ExtractComponents(vh);

    % Combine the sperm head and sperm tail information===============================================================    
    xi = [xh1;xt1q;xh2;xt2q;xh3;xt3q];
    Xi = [Xh1;Xt1q;Xh2;Xt2q;Xh3;Xt3q];
    vi = [vh1;vt1q;vh2;vt2q;vh3;vt3q];

%%  Position of the sperm head and tail in the lab frame of reference

    xx0=ApplyRotationMatrix(B,xi); % x-x0
    xs=TranslatePoints(xx0,x0);
    vs=ApplyRotationMatrix(B,vi);
    Xs=ApplyRotationMatrix(B,Xi);
    Xs=TranslatePoints(Xs,x0);
    
    NNss=NearestNeighbourMatrix(Xs,xs,blockSize);
    
    if isempty(boundary)
        xb=[];
        Xb=[];
        vb=[];
        NNbb=sparse([]);
    else
        [xb,Xb]=boundary.fn(boundary.model);
        vb=xb*0; % boundary is stationary
        NNbb=NearestNeighbourMatrix(Xb,xb,blockSize);
    end
    
    x=MergeVectorGrids(xs,xb);
    X=MergeVectorGrids(Xs,Xb);
    v=MergeVectorGrids(vs,vb);
    
    NN=MergeNNMatrices(NNss,NNbb); % assemble nearest-neighbour matrices separately for swimmer and body
 
    % now formulate fluid dynamics problem... essentially a mobility problem 

    Ns=length(xs)/3;
    Qs=length(Xs)/3;
    Nb=length(xb)/3;
    
    N=Ns+Nb;

    [AS,~]=AssembleStokesletMatrix(x,X,x,epsilon,domain,blockSize,NN);

    AU =-kron(eye(3),[ones(Ns,1);zeros(Nb,1)]) ; % component of velocity due to translational velocity of swimmer; zero velocity of boundary
    
    AF = kron(eye(3),[sum(NNss(1:Qs,1:Ns),1) zeros(1,Nb)]); % force summation - only on swimmer
        
    [x1,x2,x3]=ExtractComponents(xx0);ze=0*x1; % component of velocity due to rotation of swimmer about x0; zero velocity of boundary
    AOm=[ze -x3 x2; zeros(Nb,3); x3 ze -x1; zeros(Nb,3); -x2 x1 ze; zeros(Nb,3)];

    [x1,x2,x3]=ExtractComponents(Xs'*NNss);ze=0*x1; % moment summation
    
    AM=[ ze zeros(1,Nb) -x3 zeros(1,Nb)  x2 zeros(1,Nb); ...
         x3 zeros(1,Nb)  ze zeros(1,Nb) -x1 zeros(1,Nb); ...
        -x2 zeros(1,Nb)  x1 zeros(1,Nb)  ze zeros(1,Nb)];  

    A=[AS AU AOm; AF zeros(3,6); AM zeros(3,6)];
    
    % rhs assembly
    rhs=[v;zeros(6,1)];

    % solve and extract f, U, Omega
    sol=A\rhs;
    f=sol(1:3*N);
    U=sol(3*N+1:3*N+3);
    Om=sol(3*N+4:3*N+6);

    dz(1:3,1)=U;
    dz(4:6,1)=cross(Om,b1);
    dz(7:9,1)=cross(Om,b2);
    
    if ~isempty(varargin)
        dz(10:9+3*N,1)=f;
    end
    
    