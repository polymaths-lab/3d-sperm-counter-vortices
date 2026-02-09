% Derivative using interpolation

% derivate any function in respect to s or t
% df/ds -> if 'a' is a the arclenght space matrix
% df/dt -> take the transpost and then the derivative..

function [dFdt] = Function_derivative_of_time(Function,tRange,var)

%%
if 1
    n_q    = 30; % NUMBER OF POINTS IN THE INTERPOLATION


    % data size
    ns=size(Function,1); % points on a curve at a given time
    nt=size(Function,2); % time points

    tRange=repmat(tRange,ns,1);   
    dFdt = zeros(size(Function));


    if strcmp(var,'t') || strcmp(var,'time')
        delta =1/n_q* [tRange(:,2)-tRange(:,1),(tRange(:,3:end)-tRange(:,1:end-2))/2,tRange(:,end)-tRange(:,end-1)];
       for i_s = 1:ns
            ti_left   = tRange(i_s,:) - delta(i_s,:);
            ti_right  = tRange(i_s,:) + delta(i_s,:);
            ti=tRange(i_s,:);
            fi=Function(i_s,:);
            f_left    = spline(ti,fi,ti_left);
            f_right   = spline(ti,fi,ti_right);
            dFdt(i_s,:) = (f_right-f_left)./(2*delta(i_s,:));                      
        end

   
    end

end

%%
% if  strcmp(var,'t')
%     Function = Function';
%     a = a';
% end
% 
% np = size(Function,1); % arclenth index
% nf = size(Function,2); % arclenth index
% 
% 
% % derivate considering that the 'a' diference is contant
% if (isscalar(a))
%     dF = [(Function(2,:)-Function(1,:))         ./ (   a );...
%           (Function(3:np,:)-Function(1:np-2,:)) ./ ( 2*a );...
%           ( Function(np,:)-Function(np-1,:) )   ./ (   a )];
% else
%     % derivate considering that the 'a' varys: a = matrix
%     
%    dF = [(Function(2,:)-Function(1,:))         ./ ( a(2,:)-a(1,:)         );...
%          (Function(3:np,:)-Function(1:np-2,:)) ./ ( a(3:np,:)-a(1:np-2,:) );...
%          ( Function(np,:)-Function(np-1,:) )   ./ ( a(np,:)-a(np-1,:)     )];
% end
% 
% % taking the transpost again to back to normal form
% if strcmp(var,'t')
%     dF = dF';
% end
