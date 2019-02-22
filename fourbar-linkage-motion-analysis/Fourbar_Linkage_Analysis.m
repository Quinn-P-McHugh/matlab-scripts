%% Fourbar Linkage Position, Velocity, and Acceleration Analysis
% Uses vector loop approach to solve for positions, velocities, and accelerations of the fourbar linkage.

% Quinn McHugh
% Intro to Mechanical Design - Section 3 
% Professor Dyer
 
%% Prepare Workspace
clear; close all; clc
 
 
%% Variables for Position and Acceleration Analysis
a = 0.2;     % Crank length [m]
b = 0.5;     % Coupler length [m]
c = 0.3;     % Rocker length [m]
d = 0.55;    % Distance between pins [m]
p = 0.5/2;   % Length to the rider [m]
 
gamma = 0*pi/180;      % Angle to the rider [rad]
omega2 = (10/3)*pi;    % Angular velocity of crank [rad/s]
alpha2 = 0;            % Angular acceleration of motor [rad/s^2]
 
 
%% Variables for Force Analysis
m2 = 0.6460; m3 = 101.4560; m4 = 0.9160;        % Masses of each link [kg]
 
I2 = 0.003231; I3 = 0.03563; I4 = 0.009000;     % Moments of inertia for each link [kg*m^2]
 
fPx  = 0; fPy  = 0;         % Force applied on link 3
r3Px = 0; r3Py = 0;         % Distance of force applied from center of mass
 
T4 = 0;                     % Torque on link 4
 
 
%% Preallocate Space for Position and Acceleration Analysis
theta2 = zeros(361,1); theta3 = zeros(361,1); theta4 = zeros(361,1);    % Angles
 
xA = zeros(361,1); yA = zeros(361,1);               % Position of point A
xB = zeros(361,1); yB = zeros(361,1);               % Position of point B
xR = zeros(361,1); yR = zeros(361,1);               % Position of the rider
 
vAx = zeros(361,1); vAy = zeros(361,1);             % Velocity of point A
vRx = zeros(361,1); vRy = zeros(361,1);             % Velocity of the rider
 
omega3 = zeros(361,1); omega4 = zeros(361,1);       % Angular velocities
 
aBx = zeros(361,1); aBy = zeros(361,1);             % Acceleration of point A
aRx = zeros(361,1); aRy = zeros(361,1);             % Acceleration of the rider
 
alpha3 = zeros (361,1); alpha4= zeros(361,1);       % Angular accelerations
 
 
%% Preallocate Space for Force Analysis
a2x = zeros(361,1); a2y = zeros(361,1);         % Acceleration of crank's center of mass
a3x = zeros(361,1); a3y = zeros(361,1);         % Acceleration of coupler's center of mass
a4x = zeros(361,1); a4y = zeros(361,1);         % Acceleration of rocker's center of mass
 
fAx = zeros(361,1); fAy = zeros(361,1);         % Force on point A
fBx = zeros(361,1); fBy = zeros(361,1);         % Force on point B
fCx = zeros(361,1); fCy = zeros(361,1);         % Force on point C
fDx = zeros(361,1); fDy = zeros(361,1);         % Force on point D
 
T2 = zeros(361,1);                              % Torque on crank
 
 
%% Main Program
for i = 1:361
    
    
%% Position, Velocity, and Acceleration Analysis of Fourbar Linkage
    theta2(i) = (i-1)*pi/180;       % Crank angle
    
% Calculate variables dependent upon theta2
    r     = d - a*cos(theta2(i));
    s     = a*sin(theta2(i));
    f2    = r^2 + s^2;
    delta = acos((b^2+c^2-f2)/(2*b*c));
    g     = b - c*cos(delta);
    h     = c*sin(delta);
    
% Calculate angles of coupler and rocker
    theta3(i) = atan2(h*r - g*s, g*r + h*s);
    theta4(i) = theta3(i) + delta;
 
% Calculate X and Y position for points A, B, and the rider
    xA(i) = a*cos(theta2(i));                   % x coordinate of point A 
    yA(i) = a*sin(theta2(i));                   % y coordinate of point A 
    xB(i) = c*cos(theta4(i)) + d;               % x coordinate of point B 
    yB(i) = c*sin(theta4(i));                   % y coordinate of point B 
    
    xR(i) = xA(i) + p*cos(theta3(i)+gamma);     % x coordinate of the rider
    yR(i) = yA(i) + p*sin(theta3(i)+gamma);     % y coordinate of the rider
 
% Calculate A matrix and b vector to find velocity
    A_Mat = [b*sin(theta3(i)) -c*sin(theta4(i));
             b*cos(theta3(i)) -c*cos(theta4(i))];
 
    b_Vec = -omega2*a*[sin(theta2(i)); 
                       cos(theta2(i))];
    
% Solve for velocity vector as [omega3; omega4]
    x_Vec = A_Mat\b_Vec;
 
    omega3(i) = x_Vec(1);
    omega4(i) = x_Vec(2);
    
% Calculate X and Y velocity for point A and the rider
    vAx(i) = -omega2*a*sin(theta2(i));                      % x velocity of point A
    vAy(i) =  omega2*a*cos(theta2(i));                      % y velocity of point A
    
    vRx(i) = vAx(i) - omega3(i)*p*sin(theta3(i)+gamma);     % x velocity of the rider
    vRy(i) = vAy(i) + omega3(i)*p*cos(theta3(i)+gamma);     % y velocity of the rider
 
% Calculate C matrix and d vector to find acceleration
    C_Mat = [-b*sin(theta3(i))  c*sin(theta4(i));
              b*cos(theta3(i)) -c*cos(theta4(i))];
     
    d_Vec = [a*alpha2*sin(theta2(i))+omega2*omega2*a*cos(theta2(i))+omega3(i)*omega3(i)*b*cos(theta3(i))-omega4(i)*omega4(i)*c*cos(theta4(i)); 
            -a*alpha2*cos(theta2(i))+omega2*omega2*a*sin(theta2(i))+omega3(i)*omega3(i)*b*sin(theta3(i))-omega4(i)*omega4(i)*c*sin(theta4(i))];
      
% Solve for acceleration vector as [alpha3; alpha4]
    alpha_Vec = C_Mat\d_Vec;
 
    alpha3(i) = alpha_Vec(1);
    alpha4(i) = alpha_Vec(2);
    
% Calculate X and Y acceleration for points A and the rider
    aBx(i) = -alpha2*a*sin(theta2(i))-omega2*omega2*a*cos(theta2(i));       % x acceleration of point B
    aBy(i) =  alpha2*a*cos(theta2(i))-omega2*omega2*a*sin(theta2(i));       % y acceleration of point B
    
    aRx(i) = aBx(i) - alpha3(i)*p*sin(theta3(i)+gamma)-omega3(i)*omega3(i)*p*cos(theta3(i)+gamma);      % x acceleration of the rider
    aRy(i) = aBy(i) + alpha3(i)*p*cos(theta3(i)+gamma)-omega3(i)*omega3(i)*p*sin(theta3(i)+gamma);      % y acceleration of the rider
    
    
%% Force Analysis of Fourbar Linkage
 
% Calculate distances of point A, B, C, and D from each link's center of mass
    r2Ax = -(a/2)*cos(theta2(i));      
    r2Ay = -(a/2)*sin(theta2(i));
    r2Bx =  (a/2)*cos(theta2(i));
    r2By =  (a/2)*sin(theta2(i));
    
    r3Bx = -(b/2)*cos(theta3(i));
    r3By = -(b/2)*sin(theta3(i));
    r3Cx =  (b/2)*cos(theta3(i));
    r3Cy =  (b/2)*sin(theta3(i));
    
    r4Cx =  (c/2)*cos(theta4(i));
    r4Cy =  (c/2)*sin(theta4(i));
    r4Dx = -(c/2)*cos(theta4(i));
    r4Dy = -(c/2)*sin(theta4(i));
  
    a2x(i) = -alpha2*(a/2)*sin(theta2(i))-omega2*omega2*(a/2)*cos(theta2(i));       % x acceleration of crank's center of mass
    a2y(i) =  alpha2*(a/2)*cos(theta2(i))-omega2*omega2*(a/2)*sin(theta2(i));       % y acceleration of crank's center of mass
    
    a3x(i) = aRx(i);        % x acceleration of coupler's center of mass
    a3y(i) = aRy(i);        % y acceleration of coupler's center of mass
 
    a4x(i) = -alpha4(i)*(c/2)*sin(theta4(i))-omega4(i)*omega4(i)*(c/2)*cos(theta4(i));      % x acceleration of rocker's center of mass
    a4y(i) =  alpha4(i)*(c/2)*cos(theta4(i))-omega4(i)*omega4(i)*(c/2)*sin(theta4(i));      % y acceleration of rocker's center of mass
 
% Calculate S matrix and t vector to find force and torque
    S_Mat = [ 1     0     1     0     0     0     0    0    0;
              0     1     0     1     0     0     0    0    0;
             -r2Ay  r2Ax -r2By  r2Bx  0     0     0    0    1;
              0     0    -1     0     1     0     0    0    0;
              0     0     0    -1     0     1     0    0    0;
              0     0     r3By -r3Bx -r3Cy  r3Cx  0    0    0;
              0     0     0     0    -1     0     1    0    0;
              0     0     0     0     0    -1     0    1    0;
              0     0     0     0     r4Cy -r4Cx -r4Dy r4Dx 0];
          
    t_Vec = [m2*a2x(i);
             m2*a2y(i);
             I2*alpha2;
             m3*a3x(i)-fPx;
             m3*a3y(i)-fPy;
            -r3Px*fPy+r3Py*fPx+I3*alpha3(i);
             m4*a4x(i);
             m4*a4y(i);
             I4*alpha4(i)-T4];
      
% Solve for force and torque vector
    f_Vec = S_Mat\t_Vec;
 
    fAx(i) = f_Vec(1);
    fAy(i) = f_Vec(2);
    fBx(i) = f_Vec(3);
    fBy(i) = f_Vec(4);
    fCx(i) = f_Vec(5);
    fCy(i) = f_Vec(6);
    fDx(i) = f_Vec(7);
    fDy(i) = f_Vec(8);
    T2(i)  = f_Vec(9);
    
    
end
 
 
%% Plot Position and Acceleration Analysis
ang = 60; % This is the crank angle where we plot a "snapshot" 
figure; hold on; grid on; axis equal
patch([xA(ang) xB(ang) xR(ang)],[yA(ang) yB(ang) yR(ang)],[0.9 1 1]) 
plot(xR,yR,'LineWidth',3,'Color','m')
xlabel('X Position [m]'); ylabel('Y Position [m]')
title('Trajectory of Rider for One Rotation of Crank')
plot([0 xA(ang) xB(ang) d],...
     [0 yA(ang) yB(ang) 0],'LineWidth',2,'Color','k')
plot(0,0,'o','Markersize',10,'MarkerFaceColor','k','Color','k')
plot(d,0,'o','Markersize',10,'MarkerFaceColor','k','Color','k')
plot(xR(ang),yR(ang),'o','Markersize',10,'MarkerFaceColor','r','Color','r')
 
figure; hold on; axis([0 360 -4 3]); grid on
set(gca,'XTick',0:60:360)
plot(theta2*180/pi,aRx*0.101971621,'LineWidth',2,'Color','b')
plot(theta2*180/pi,aRy*0.101971621,'LineWidth',2,'Color','m')
legend('X Acceleration','Y Acceleration','Location','Northeast')
title('X and Y Acceleration of Rider for One Rotation of Crank')
xlabel('Crank Angle [deg]'); ylabel('Acceleration [g]')
 
 
%% Plot Force Analysis
figure; hold on; axis([0 360 -450 250]); grid on
set(gca,'XTick',0:60:360)
plot(theta2*180/pi,T2,'LineWidth',2,'Color','b')
title('Torque Required to Spin the Crank at 100RPM v.s. Crank Angle')
xlabel('Crank Angle [deg]'); ylabel('Torque [N*m]')

