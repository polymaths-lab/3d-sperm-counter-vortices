function outputs=DKActSpermWave_modified(xii,t,outopt,alpha_temp,varargin)

% generates 'activated'waveform used by Dresdner & Katz 1981

% if varargin is non-empty, wavenumber is read in from first entry

% ds is arclength metric

%a=0.1087;
%b=0.0543;

if ~isempty(varargin)
    k=varargin{1}.k;
    phase=varargin{1}.phase;
else
    k=2*pi;
    phase=0;
end

t=t+phase;

% set up in 'wave frame' first
%yii=(a*xii+b).*sin(k*xii-t)-b*sin(-t);

% rotation angle into body frame (so flagellum is tangential to head)
%th=atan(a*sin(-t)+k*b*cos(-t));

% rotate into body frame
%xi= cos(th)*xii+sin(th)*yii;
%yi=-sin(th)*xii+cos(th)*yii;

a=0.2;
w=1;
B=0;
arf=alpha_temp;  %0; 
n=1.0;
xi= xii;
% For pusher
yi=a*(xi.^n).*(cos(k*xi-w*t)+B); 
zi=-arf*a*(xi.^n).*sin(k*xi-w*t); 
% For puller
%yi=a*(xi.^n).*(cos(k*xi+w*t)+B); 
%zi=-arf*a*(xi.^n).*sin(k*xi+w*t); 

switch outopt
    case 'xy'
        outputs=[xi,yi,zi];
    case 'ds'
        % For pusher
        outputs=sqrt(1+(a*n*(xi.^(n-1)).*(cos(k*xi-w*t)+B)-k*a*(xi.^n).*sin(k*xi-w*t)).^2+(arf*a*n*(xi.^(n-1)).*sin(k*xi-w*t)+k*arf*a*(xi.^n).*cos(k*xi-w*t)).^2); 
        % For puller
        %outputs=sqrt(1+(a*n*(xi.^(n-1)).*(cos(k*xi+w*t)+B)-k*a*(xi.^n).*sin(k*xi+w*t)).^2+(arf*a*n*(xi.^(n-1)).*sin(k*xi+w*t)+k*arf*a*(xi.^n).*cos(k*xi+w*t)).^2);    
end