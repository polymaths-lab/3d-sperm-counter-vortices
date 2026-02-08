function hp=PlotMultiSpermSolution(fig,swimmer,plotSwitch,t,z,m,varargin)

if ~isempty(varargin)
    mksz=varargin{1};
else
    mksz=1;
end

Nsw=length(swimmer);

if Nsw==1
    swimmertemp=swimmer;
    clear swimmer;
    swimmer{1}=swimmertemp;
end

x0 = cell(1,Nsw);
B = cell(1,Nsw); b1 = cell(1,Nsw); b2 = cell(1,Nsw); b3 = cell(1,Nsw);
xi = cell(1,Nsw); vi = cell(1,Nsw); Xi = cell(1,Nsw);
x = cell(1,Nsw); x1 = cell(1,Nsw); x2 = cell(1,Nsw); x3 = cell(1,Nsw);
X = cell(1,Nsw); X1 = cell(1,Nsw); X2 = cell(1,Nsw); X3 = cell(1,Nsw);
O = cell(1,Nsw);

for n=1:Nsw
    x0{n}=z(m,       n:Nsw:n+2*Nsw);             %3*(n-1)+1:3*n);
    b1{n}=z(m, 3*Nsw+n:Nsw:n+5*Nsw);       %:3*Nsw+3*n);
    b2{n}=z(m, 6*Nsw+n:Nsw:n+8*Nsw);       %+3*(n-1)+1:6*Nsw+3*n);
    b3{n}=cross(b1{n},b2{n});
    B{n}=[b1{n}(:) b2{n}(:) b3{n}(:)];
    O{n}=x0{n}; %[x0{n}(1) x0{n}(2) x0{n}(3)];
    
    
    [xi{n},vi{n},Xi{n}]=swimmer{n}.fn(t(m),swimmer{n}.model);
    x{n}=ApplyRotationMatrix(B{n},xi{n});
    X{n}=ApplyRotationMatrix(B{n},Xi{n});
    x{n}=TranslatePoints(x{n},O{n});
    X{n}=TranslatePoints(X{n},O{n});
    [x1{n},x2{n},x3{n}]=ExtractComponents(x{n});
    [X1{n},X2{n},X3{n}]=ExtractComponents(X{n});
%----------------------
   % M=swimmer{n}.model.nh;%6*swimmer{n}.model.nh*swimmer{n}.model.nh;
    %kT=length(X1{n})-M;
    k=length(X1{n})-50;%round(kT/2);
    %round(length(X1t{n})/2);%定位 标记点p
%--------------------------
end

figure(fig);hold on;
for n=1:Nsw
    switch plotSwitch
        case 'f' %force points
            hp=plot3(x1{n},x2{n},x3{n},'b.');set(hp,'markersize',mksz);
            %-----------------------------------------------------------
        case 'q' %quadrature points
            hp=plot3(X1{n},X2{n},X3{n},'r.');set(hp,'markersize',mksz);
        case 'proXY' %quadrature points投影
            [j k]=size(X3{n});
            Z{n}=zeros(j,k);
            hp=plot3(X1{n},X2{n},Z{n},'r.');set(hp,'markersize',mksz);
        case 'proXZ' %quadrature points投影
            [j k]=size(X2{n});
            Y{n}=ones(j,k)*1;
            hp=plot3(X1{n},Y{n},X3{n},'r.');set(hp,'markersize',mksz);
        case 'proYZ' %quadrature points投影
            [j k]=size(X1{n});
            X{n}=ones(j,k)*1;
            hp=plot3(X{n},X2{n},X3{n},'r.');set(hp,'markersize',mksz);
            %-------------------------------------------------------------
        case 'p' %标记鞭毛中间点（方便画出该点的轨迹）
            hp=[X1{n}(k),X2{n}(k),X3{n}(k)];%修改
        case 'pXY' %标记鞭毛中间点（方便画出该点的轨迹）投影
            hp=[X1{n}(k),X2{n}(k),0];
        case 'pXZ' %标记鞭毛中间点（方便画出该点的轨迹）投影
            hp=[X1{n}(k),1,X3{n}(k)];
        case 'pYZ' %标记鞭毛中间点（方便画出该点的轨迹）投影
            hp=[1,X2{n}(k),X3{n}(k)]
           %-------------------------------------------------------------
        case 'a' %force and quadrature points
            hp=plot3(x1{n},x2{n},x3{n},'b.');set(hp,'markersize',mksz);
            hp=plot3(X1{n},X2{n},X3{n},'r.');set(hp,'markersize',mksz);
    end
end
