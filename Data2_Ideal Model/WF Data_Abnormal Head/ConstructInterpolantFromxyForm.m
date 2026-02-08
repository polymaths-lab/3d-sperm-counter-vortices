function F=ConstructInterpolantFromxyForm(s,t,xyWaveFn,varargin)

% constructs a griddedInterpolant for sperm waveforms of xy form
%
% (s,t) are produced by ndgrid
% xyWaveFn - function to give wave in xy form (as D&K)
%
x = zeros(size(s,1),size(s,2));
y = x;
z = x; %修改
for k=1:size(s,2)
    for j=1:size(s,1)
        if isempty(varargin)
            [x(j,k),y(j,k),z(j,k)]=FindxyAtArclength(s(j,k),t(j,k),xyWaveFn);%修改 
        else
            [x(j,k),y(j,k),z(j,k)]=FindxyAtArclength(s(j,k),t(j,k),xyWaveFn,...
                varargin{1});%修改
        end
    end
end

F{1}=griddedInterpolant(s,t,x,'spline');
F{2}=griddedInterpolant(s,t,y,'spline');
F{3}=griddedInterpolant(s,t,z,'spline');%修改