function [x,y,z]=CalcxyFromPlanarInterp(s,t,F)

% calculates the x and y coordinates from arclength s (array)
% time t (scalar), gridded-interpolants F{1}, F{2}

x=F{1}(s,t);
y=F{2}(s,t);
z=F{3}(s,t);


