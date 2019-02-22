% Triple_Pendulum_Plot.m
% makes primitive animation of triple pendulum

function Triple_Pendulum_Plot(pend,q0,t,q)

L1 = pend.L1; L2 = pend.L2; L3 = pend.L3;

% nStep = length(t);  % Total number of time steps
% nShot = 200;        % Number of "snapshots" to plot
% nStrobe = floor(nStep/nShot);   % Number of time steps between plots
% 
% % Open figure window
% pend_plot = figure;
% set(pend_plot,'Name','Triple Pendulum Movie',...
%     'Units','Normalized','OuterPosition',[0.0 0.5 0.5 0.5])
% 
% % Create initial plot in initial configuration
% OR = [0; 0;];
% M1 = [L1*cos(q0(1));                                    L1*sin(q0(1))]/0.0254;
% M2 = [L1*cos(q0(1)) + L2*cos(q0(2));                    L1*sin(q0(1)) + L2*sin(q0(2))]/0.0254;
% M3 = [L1*cos(q0(1)) + L2*cos(q0(2)) + L3*cos(q0(3));    L1*sin(q0(1)) + L2*sin(q0(2)) + L3*sin(q0(3))]/0.0254;
% 
% Pend1 = [OR M1];    % Upper pendulum
% Pend2 = [M1 M2];    % Middle pendulum
% Pend3 = [M2 M3];    % Lower pendulum
% 
% hold on; grid on; axis equal;
% 
% for i = 1:nStrobe:nStep
%     M1 = [L1*cos(q(i,1));                        L1*sin(q(i,1))]/0.0254; % position of M1
%     M2 = [L1*cos(q(i,1)) + L2*cos(q(i,2));...
%           L1*sin(q(i,1)) + L2*sin(q(i,2))]/0.0254; % position of M2
%     M3 = [L1*cos(q(i,1)) + L2*cos(q(i,2)) + L3*cos(q(i,3));...
%           L1*sin(q(i,1)) + L3*sin(q(i,2)) + L3*sin(q(i,3))]/0.0254; % position of M3
% 
%     Pend1 = [OR M1]; % upper pendulum
%     Pend2 = [M1 M2]; % middle pendulum
%     Pend3 = [M2 M3]; % lower pendulum
% 
%     plot(Pend1(1,:), Pend1(2,:), 'LineWidth',2,'Color','y')
%     plot(Pend2(1,:), Pend2(2,:), 'LineWidth',2,'Color','b')
%     plot(Pend3(1,:), Pend3(2,:), 'LineWidth',2,'Color','g')
% 
%     plot(M1(1), M1(2),'o','MarkerSize',3,'MarkerFaceColor','b')
%     plot(M2(1), M2(2),'o','MarkerSize',5,'MarkerFaceColor','k')
%     plot(M3(1), M3(2),'o','MarkerSize',5,'MarkerFaceColor','r')
% 
%     pause(0.02)
% end

% Trajectory of mass 1
figure
plot((L1*cos(q(:,1)))/0.0254, (L1*sin(q(:,1)))/0.0254,'color','m')
title('Trajectory of First Pendulum')
xlabel('x displacement (in)')
ylabel('y displacement (in)')
axis equal; grid on; 

% Trajectory of mass 2
figure 
plot((L1*cos(q(:,1)) + L2*cos(q(:,2)))/0.0254, (L1*sin(q(:,1)) + L2*sin(q(:,2)))/0.0254,'color','c')
title('Trajectory of Second Pendulum')
xlabel('x displacement (in)')
ylabel('y displacement (in)')
axis equal; grid on;

% Trajectory of mass 3
figure
plot((L1*cos(q(:,1)) + L2*cos(q(:,2)) + L3*cos(q(:,3)))/0.0254, (L1*sin(q(:,1)) + L2*sin(q(:,2)) + L3*sin(q(:,3)))/0.0254,'color','r')
title('Trajectory of Third Pendulum')
xlabel('x displacement (in)')
ylabel('y displacement (in)')
axis equal; grid on;


