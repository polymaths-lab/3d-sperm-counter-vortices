function outputs=DKActSpermWave(xii,t,outopt,varargin)

% generates 3d waveform modified based on the model in 'Smith D J, 2009, Human sperm accumulation',
% which originates from the 2d waveform model by 'Dresdner & Katz, 1981,Relationships of Mammalian'.

% if varargin is non-empty, wavenumber is read in from first entry

% ds is arclength metric

a=0.2;
if ~isempty(varargin)
    k=varargin{1}.k;
    phase=varargin{1}.phase;
else
    k=2*pi;
    phase=0;
end
t=t+phase;
w=1;
B=0;
arf= 0.6; 
n=1.0;

xi= xii;
yi=a*(xi.^n).*(cos(k*xi-w*t)+B); 
zi=-arf*a*(xi.^n).*sin(k*xi-w*t); 


switch outopt
    case 'xy'
        outputs=[xi,yi,zi];
    case 'ds'
        outputs=sqrt(1+(a*n*(xi.^(n-1)).*(cos(k*xi-w*t)+B)-k*a*(xi.^n).*sin(k*xi-w*t)).^2+(arf*a*n*(xi.^(n-1)).*sin(k*xi-w*t)+k*arf*a*(xi.^n).*cos(k*xi-w*t)).^2);            
end