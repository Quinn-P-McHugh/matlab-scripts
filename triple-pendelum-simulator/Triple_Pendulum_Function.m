% Triple_Pendulum_Function.m
% Gives differential equation for triple pendulum

function df = Triple_Pendulum_Function(~,q,pend)

% Dimensions
L1 = pend.L1; L2 = pend.L2; L3 = pend.L3;
m1 = pend.m1; m2 = pend.m2; m3 = pend.m3;
g = pend.g;

% Calculated inertias
I123 = (m2 + m3)*L1*L2;
I13 = m3*L1*L3;
I23 = m3*L2*L3;
J1 = (m1 + m2 + m3)*L1;
J2 = (m2 + m3)*L2;
J3 =  m3*L3;

df(1,1) = q(4);
df(2,1) = q(5);
df(3,1) = q(6);
df(4,1) = -I123*(q(5)^2)*sin(q(1)-q(2)) - I13*(q(6)^2)*sin(q(1)-q(3)) - J1*g*cos(q(1));
df(5,1) = -I123*(q(4)^2)*sin(q(2)-q(1)) - I23*(q(6)^2)*sin(q(2)-q(3)) - J2*g*cos(q(2));
df(6,1) = -I13*(q(4)^2)*sin(q(3)-q(1)) - I23*(q(5)^2)*sin(q(3)-q(2)) - J3*g*cos(q(3));