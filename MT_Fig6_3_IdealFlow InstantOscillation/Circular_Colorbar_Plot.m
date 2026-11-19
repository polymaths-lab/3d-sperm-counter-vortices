function Circular_Colorbar_Plot

clc
%% Plot circular colorbar. -- reference code see 'report35' (from Hermes).
load idx_flow.mat;


N = length(idx_flow);
%N = 256;                                % Number Of Colours In Colour Wheel
th = linspace(pi/2, -3*pi/2, N);              % Angles
xy = @(a,r) [r.*cos(a); r.*sin(a)];     % Function: Calculate (x,y) Vectors
r_ext = 0.7;                            % Outer Radius
r_int = 0.3;                            % Inner Radius
C(1,:,:) = xy(th,r_ext);                % Outer Circle
C(2,:,:) = xy(th,r_int);                % Inner Circle
C(3,:,:) = xy(th,1);
%c = colormap(parula(N));                   % Set ?colormap?'
c_temp = colormap(jet(70)); 
c = [c_temp([1 16 30],:); c_temp([70 60 50],:)];

figure(1)
set(gcf,'color','w') %
clf;hold on;
grid on
axis square off

C1 = squeeze(C(1,:,:));                 % Reduce ?C1? Dimensions
plot(C1(1,:), C1(2,:))
C2 = squeeze(C(2,:,:));                 % Reduce ?C2? Dimensions
plot(C2(1,:), C2(2,:))
C3 = squeeze(C(3,:,:));
for k1 = 1:size(C2,2)
    col_idx = idx_flow(k1);
    %if col_idx==1 || col_idx==2 || col_idx==3
    plot([C1(1,k1), C2(1,k1)], [C1(2,k1), C2(2,k1)], 'Color',c(col_idx,:), 'LineWidth',0.25)
    k2 = max([mod(k1+1,size(C2,2)), 1]);
    patch([C1(1,k1), C1(1,k2), C2(1,k2), C2(1,k1)], [C1(2,k1), C1(2,k2), C2(2,k2), C2(2,k1)], c(col_idx,:), 'EdgeColor','none') 
    if mod(k1,100)==0
        plot([C2(1,k1), C3(1,k1)], [C2(2,k1), C3(2,k1)], 'k:', 'LineWidth',10)
    end
    %end
end
axis([-1,1,-1,1]*1.11)


return
fontsize1 = 160;

annotation('textbox',...
    [0.66 0.5 0.0525714285714286 0.0690476190476196],...
    'String','$0$',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',fontsize1,...
    'FitBoxToText','off');

annotation('textbox',...
    [0.5 0.87 0.0999999999999993 0.001],...
    'String','$\pi/2$',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',fontsize1,...
    'FitBoxToText','off');

% Create textbox
annotation('textbox',...
    [0.35 0.58 0.0999999999999993 0.001],...
    'String','$\pi$',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',fontsize1,...
    'FitBoxToText','off');

annotation('textbox',...
[0.49 0.24 0.0999999999999994 0.000999999999999945],...
    'String','$3\pi/2$',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',fontsize1,...
    'FitBoxToText','off');


% Create arrow
annotation('arrow',[0.7 0.65],...
    [0.678571428571428 0.773109243697478],'LineWidth',5,'HeadStyle','cback2','HeadSize',30);
