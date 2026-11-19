function [dFdt] = Function_derivative_of_time(Function,tRange)


ns=size(Function,1); 
%nt=size(Function,2); % time points

tRange=repmat(tRange,ns,1);   
dFdt = zeros(size(Function));
n_q = 30; 
delta =1/n_q* [tRange(:,2)-tRange(:,1),(tRange(:,3:end)-tRange(:,1:end-2))/2,tRange(:,end)-tRange(:,end-1)];

for i_s = 1:ns
    ti_left = tRange(i_s,:) - delta(i_s,:);
    ti_right = tRange(i_s,:) + delta(i_s,:);
    ti = tRange(i_s,:);
    fi = Function(i_s,:);
    f_left = spline(ti,fi,ti_left);
    f_right = spline(ti,fi,ti_right);
    dFdt(i_s,:) = (f_right-f_left)./(2*delta(i_s,:));                      
end

   




