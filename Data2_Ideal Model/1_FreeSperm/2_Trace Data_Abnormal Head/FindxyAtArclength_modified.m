function [x,y,z]=FindxyAtArclength_modified(s,t,xyWaveFn,alpha_temp,varargin)

% for waveforms not prescribed in terms of arclength, this function
% takes arclength coordinate input and outputs x,y coordinates at that
% arclength

if isempty(varargin) 
    Resi=@(xii) s-integral(@(xiii) xyWaveFn(xiii,t,'ds',alpha_temp),0,xii);
else
    Resi=@(xii) s-integral(@(xiii) xyWaveFn(xiii,t,'ds',alpha_temp,varargin{1}),0,xii);
end

options=optimoptions(@fsolve,'Display','off');

xii0=s;


xii=fsolve(Resi,xii0,options);

if isempty(varargin)
    outputs=DKActSpermWave_modified(xii,t,'xy',alpha_temp);
else
    outputs=DKActSpermWave_modified(xii,t,'xy',alpha_temp,varargin{1});
end

x=outputs(:,1);
y=outputs(:,2);
z=outputs(:,3);
