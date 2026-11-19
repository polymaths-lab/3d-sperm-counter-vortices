function dw = ODE_SolveFrenetEquations(curvature,torsion,w,s,sspan)

%% solve non-stiff differential equations (built-in function: ode45) of the Frenet frame

delta_s = (sspan(end)-sspan(1))/(length(sspan)-1);
ind_s = find( abs(sspan-s) < delta_s );
ind_s = ind_s (1);

h = [curvature(ind_s), torsion(ind_s)];
F = diag(h,1) - diag(h,-1);
A = kron(F, eye(3));

dw = A*w;