% Triple_Pendulum_Mass
% Solves for mass matrix of triple pendulum

function M = Triple_Pendulum_Mass(~,q,pend)

% Dimensions
L1 = pend.L1; L2 = pend.L2; L3 = pend.L3;
m1 = pend.m1; m2 = pend.m2; m3 = pend.m3;

% Calculated inertias
I1 =   (m1 + m2 + m3)*(L1^2);
I2 =   (m2 + m3)*(L2^2);
I3 =    m3*(L3^2);
I123 = (m2 + m3)*L1*L2;
I13 =   m3*L1*L3;
I23 =   m3*L2*L3;

M = [1 0 0 0                     0                     0
     0 1 0 0                     0                     0
     0 0 1 0                     0                     0
     0 0 0 I1                    I123*cos(q(1)-q(2))   I13*cos(q(1)-q(3));
     0 0 0 I123*cos(q(2)-q(1))   I2                    I23*cos(q(2)-q(3));
     0 0 0 I13*cos(q(3)-q(1))    I23*cos(q(3)-q(2))    I3];