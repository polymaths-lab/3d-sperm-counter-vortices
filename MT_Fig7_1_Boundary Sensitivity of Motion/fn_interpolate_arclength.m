function   [xi,yi,zi,ai] = fn_interpolate_arclength(x,y,z,n_points,varargin)

% n_points = np number of points desired
%np   = size(xi(srange,:),1);
nt   = size(x,2);
%Interpolate arclength again

if ~isempty(varargin)
    pp=varargin{1};
else
    pp = 0.9999;% 1= interpolation without cubic spline
end

trange = [1:nt];
%xi = zeros(n_points,trange(end)-trange(1)+1);
xi = zeros(n_points,nt);
yi = xi;
zi = xi;


for j = trange
    clear xtemp ytemp ztemp s_temp
      xtemp = x(:,j);
      ytemp = y(:,j);
      ztemp = z(:,j);
      
      sold = arclength(xtemp, ytemp, ztemp);
      snew = linspace(0,sold(end),n_points)';

xi(:,j-trange(1)+1) = fnval(csaps(sold,xtemp,pp),snew);
yi(:,j-trange(1)+1) = fnval(csaps(sold,ytemp,pp),snew);
zi(:,j-trange(1)+1) = fnval(csaps(sold,ztemp,pp),snew);

end

% EVALUATE NEW ARCLENGTH
ai = arclength(xi,yi,zi);