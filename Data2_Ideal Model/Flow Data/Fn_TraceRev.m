function np_rev_1 = Fn_TraceRev(x0,y0,z0,np_beat)
% This function measures the time point numbers per trace ([x,y,z]) revolution, for idealized sperm models with 0<alpha<1.

    nt0 = length(x0);

    % nspiral: number of spirals/ revolutions. 
    nspiral_0 = 5; % This is a casual value, with which we can obtain a preliminary aligned trajectory, and further calculate the revolution information.
    [np_rev_0, ~]= fn_rev_TrajType_HR_TR_SS_HL(x0,y0,z0,np_beat,nspiral_0); % This function originates from rxm's paper on AS, 2025.
    nt1 = min(round(3*np_rev_0),nt0);
    x1 = x0(1:nt1);
    y1 = y0(1:nt1);
    z1 = z0(1:nt1);

    
    nspiral_1 = nt1/np_rev_0;
    [np_rev_1, ~]= fn_rev_TrajType_HR_TR_SS_HL(x1,y1,z1,np_beat,nspiral_1);
    
    

   
    

    