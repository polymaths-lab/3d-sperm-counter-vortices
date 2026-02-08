function [np_rev, pit]= fn_rev_TrajType_HR_TR_SS_HL(x0,y0,z0,np_beat,nspiral)



[~, Traj_X_new_1] = Align_traj_with_AP_CA(x0,y0,z0,nspiral);
x0_new=Traj_X_new_1(:,1);
y0_new=Traj_X_new_1(:,2);
z0_new=Traj_X_new_1(:,3); 



%% Step 1: find local maximal points, by connecting which we can get revolution cycles.

d=(y0_new.^2+z0_new.^2).^0.5; 
for i=2:length(d)-1
    if d(i)>d(i-1) && d(i)>d(i+1)
        Ind_temp = i;
        break
    end
end
count_wig_max = 1;
Ind_wig_max(count_wig_max) = Ind_temp;

for i_new = Ind_temp:np_beat:length(d)
    count_wig_max = count_wig_max + 1;
    Ind_wig_max(count_wig_max) = i_new;
end
x0_new_wig_max = x0_new(Ind_wig_max); 
y0_new_wig_max = y0_new(Ind_wig_max); 
z0_new_wig_max = z0_new(Ind_wig_max);


%% Step 2: calculate unwrapped azimuth angle the trajectory goes through, every 2pi representing 1 revolution cycle.

azi_rev_temp0 = cart2pol(y0_new_wig_max,z0_new_wig_max); 
azi_rev_temp1 = nan(length(azi_rev_temp0)-1,1);  
for i_azi = 1:length(azi_rev_temp0)-1
    if sign(azi_rev_temp0(i_azi+1))*sign(azi_rev_temp0(i_azi))<0
        delta_1 = abs(azi_rev_temp0(i_azi+1)-azi_rev_temp0(i_azi));
        delta_2 = 2*pi-delta_1;
        azi_rev_temp1(i_azi) = min(delta_1,delta_2);
    else
        azi_rev_temp1(i_azi) = abs(azi_rev_temp0(i_azi+1)-azi_rev_temp0(i_azi));
    end
end

% unwrapped azimuth
azi_rev_temp2 = cumsum(azi_rev_temp1); 
azi_rev_temp3 =nan(length(azi_rev_temp2)+1,1); 
azi_rev_temp3(1) = 0;
azi_rev_temp3(2:end) = azi_rev_temp2; 



%% Step 3: get revolution information.

nrev_temp = min(3,azi_rev_temp3(end)/(2*pi));
ind_new_path(1) = Ind_wig_max(1);
nrev_check= 0;
for i = 1:length(azi_rev_temp3)    
    if azi_rev_temp3(i)>=2*pi 
        nrev_check = nrev_check+1;
        ind_new_path(nrev_check+1) = Ind_wig_max(i);          
        if nrev_check >= floor(nrev_temp)
            break
        end  
        azi_rev_temp3 = azi_rev_temp3-2*pi; 
    end
end
  
    
np_rev = mean(diff(ind_new_path));
pit = mean(abs(diff(x0_new(ind_new_path))));






