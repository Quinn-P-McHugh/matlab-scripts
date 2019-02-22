% Triple_Pendulum_Solver.m
% Use MATLAB's RK diff eq solver to solve
% triple-pendulum problem

clear variables; close all; clc;

pend.g = 9.81;    % Gravitational constant (m/s2)
pend.m1 = 2.8;    % Mass of pendulum 1 (kg)
pend.m2 = 4.1;    % Mass of pendulum 2 (kg)
pend.m3 = 3.9;    % Mass of pendulum 3 (kg)
pend.L1 = 15.4*0.0254;      % Length of pendulum 1 (m)
pend.L2 = 20.9*0.0254;      % Length of pendulum 2 (m)
pend.L3 = 22.0*0.0254;      % Length of pendulum 3 (m)

q0 = [28.4*(pi/180) -23.9*(pi/180) -2.5*(pi/180) 0 0 0];  % Initial conditions (first 3 are displacement, second 3 are velocity)

dT = [0 5];     % Simulation time (init, final)

tol = 1e-4 * [1 1 1 1 1 1];     % Error tolerance for solver

options = odeset('RelTol',1e-4,'AbsTol',tol,...
                 'Mass',@(t,q) Triple_Pendulum_Mass(t,q,pend));
             
[T,Q] = ode45(@(t,q) Triple_Pendulum_Function(t,q,pend), dT, q0, options);

Triple_Pendulum_Plot(pend, q0, T, Q);