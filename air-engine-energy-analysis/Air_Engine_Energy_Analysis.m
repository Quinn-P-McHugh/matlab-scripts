%% Air Engine Energy Analysis
% Quinn McHugh, Carly Jorgenson, Alex Lindner, Leslie Maier, Robert
% Livingston
% Dynamics - Section 2
% Professor Osta

% Uses RPM and time data to find angular velocity, angular acceleration, 
% and energy losses of an air engine's flywheel. by dropping a weight 
% from a height of approximately 16.5ft. The weight is attached to the air 
% engine's flywheel using string.


%% User Variables
m_Weight = 0.5721; % [kg] Mass of one weight
m_Holder = 0.4272; % [kg] Mass of weight holder
m_Carabiner = 0.0204; % [kg] Mass of cararbiner (Attaches string to weight 
% holder)

g = 9.80665; % [m/s^2] Acceleration of gravity

r_BS = 0.039497/2; % [m] Radius of the big spool
r_LS = 0.019177/2; % [m] Radius of the little spool

h = 5.0292; % [m] Height of the weight when dropped


%% Finding Angular Velocity Given RPM and Time

% Calculates RPM and time for the big spool - # of weights = 1:4
BS_RPM = [BS_1W_RPM, BS_2W_RPM, BS_3W_RPM, BS_4W_RPM]; % [rev/min]
BS_Time_ms = [BS_1W_Time, BS_2W_Time, BS_3W_Time, BS_4W_Time]; % [ms]
BS_Time = BS_Time_ms/1000; % [s]

% Calculates RPM and time for the little spool - # of weights = 2:5
LS_RPM = [LS_2W_RPM, LS_3W_RPM, LS_4W_RPM, LS_5W_RPM]; % [rev/min]
LS_Time_ms = [LS_2W_Time, LS_3W_Time, LS_4W_Time, LS_5W_Time]; % [ms]
LS_Time = LS_Time_ms/1000; % [s]

% Calculates angular velocities of big and little spool
BS_Omega = ((2*pi)/60)*BS_RPM; % [rad/s]
LS_Omega = ((2*pi)/60)*LS_RPM; % [rad/s]


%% Finding Average Angular Acceleration Using Angular Velocity

% Calculates min angular velocities for big and little spools
BS_Omega_Min = BS_Omega(30); % [rad/s]
LS_Omega_Min = LS_Omega(30); % [rad/s]

% Calculates max angular velocities for the big and little spools
BS_Omega_Max(1) = BS_Omega(143,1); % [rad/s]
BS_Omega_Max(2) = BS_Omega(98,2); % [rad/s]
BS_Omega_Max(3) = BS_Omega(93,3); % [rad/s]
BS_Omega_Max(4) = BS_Omega(76,4); % [rad/s]

LS_Omega_Max(1) = LS_Omega(161,1); % [rad/s]
LS_Omega_Max(2) = LS_Omega(117,2); % [rad/s]
LS_Omega_Max(3) = LS_Omega(114,3); % [rad/s]
LS_Omega_Max(4) = LS_Omega(102,4); % [rad/s]

% Calculates average angular accelerations of the big spool (BS_Alpha)
BS_Omega_Max_I = [143, 98, 93, 76]; % Indice (row) of each max angular velocity 
% value
BS_Time_Max = BS_Time(BS_Omega_Max_I);
BS_Alpha = (BS_Omega_Max - BS_Omega_Min) ./ (BS_Time_Max - BS_Time(30)) % [rad/s^2] The 
% difference in the max and min angular velocity values divided by the 
% corresponding time at the max angular velocity value

% Calculates average angular accelerations of the little spool (LS_Alpha)
LS_Omega_Max_I = [161, 117, 114, 102]; % Indice (row) of each max angular velocity 
% value
LS_Time_Max = LS_Time(LS_Omega_Max_I);
LS_Alpha = (LS_Omega_Max - LS_Omega_Min) ./ (LS_Time_Max - LS_Time(30)) % [rad/s^2] The 
% difference in the max and min angular velocity values divided by the 
% corresponding time at the max angular velocity value


%% Finding Moment of Inertia of Flywheel + Shaft Assembly
n = 2; % Since the little spool starts with 2 weights as opposed to 1, n 
% has an initial value of 2 and increments by 1 in the loop below

for i=1:4
    BS_MomentInertia(i) = ((i*m_Weight + m_Holder + m_Carabiner)*g*r_BS) ./ BS_Alpha(i) % [kg*m^2] Moment of Inertia of the big spool
    LS_MomentInertia(i) = ((n*m_Weight + m_Holder + m_Carabiner)*g*r_LS) ./ LS_Alpha(i) % [kg*m^2] Moment of Inertia of the little spool
    

%% Using Work-Energy Equation to Find Energy Losses

    BS_Ug(i) = (i*m_Weight + m_Holder + m_Carabiner)*g*h; % [J] Gravitational potential energy
    BS_Ke(i) = (1/2)*BS_MomentInertia(i)*(BS_Omega_Max(i) .^ 2); % [J] Kinetic energy
    BS_Kr(i) = (1/2)*(i*m_Weight + m_Holder + m_Carabiner)*((BS_Omega_Max(i)*r_BS)^2); % [J] Rotational kinetic energy
    
    LS_Ug(i) = (n*m_Weight + m_Holder + m_Carabiner)*g*h; % [J] Gravitational potential energy
    LS_Ke(i) = (1/2)*LS_MomentInertia(i)*(LS_Omega_Max(i) .^ 2); % [J] Kinetic energy
    LS_Kr(i) = (1/2)*(n*m_Weight + m_Holder + m_Carabiner)*((LS_Omega_Max(i)*r_LS)^2); % [J] Rotational kinetic energy
    
    BS_Losses(i) = -BS_Ug(i) + BS_Ke(i) + BS_Kr(i) % [J] Energy losses of the big spool system
    LS_Losses(i) = -LS_Ug(i) + LS_Ke(i) + LS_Kr(i) % [J] Energy losses of the little spool system

    n = n+1; % Increments number of weights on the little spool (n) by 1
end
    
%% Standard Deviation of Moment of Inertia

LS_MomentInertia_STD = std(LS_MomentInertia)
BS_MomentInertia_STD = std(BS_MomentInertia)


%% Plot Figures

% Big Spool - Angular Velocity vs Time
figure
plot(BS_Time_ms, BS_Omega, 'LineWidth',2)
xlabel('Time [ms]','fontweight','bold','fontsize',12)
ylabel('Angular Velocity [rad/s]','fontweight','bold','fontsize',12)
for i = 1:4
    if i == 1
        BS_Legend{i} = [num2str(i),' Weight'];
    else
        BS_Legend{i} = [num2str(i),' Weights'];
    end
end
legend(BS_Legend)
        
% Little Spool - Angular Velocity vs Time
figure
plot(LS_Time_ms, LS_Omega, 'LineWidth',2)
xlabel('Time [ms]','fontweight','bold','fontsize',12)
ylabel('Angular Velocity [rad/s]','fontweight','bold','fontsize',12)
n = 2;
for i = 1:4
    LS_Legend{i} = [num2str(n),' Weights'];
    n = n + 1;
end
legend(LS_Legend)


% Angular Velocity vs Time - All Trials (Trendline must be added via Basic
% Plotting Tool)
colors = get(gca,'colororder');

figure
hold on
for i = 1:4
    % Big Spool
    plot(BS_Time_ms(20:BS_Omega_Max_I(i),i), BS_Omega(20:BS_Omega_Max_I(i),i), 'LineWidth',2,'Color',[colors(i,1),colors(i,2),colors(i,3)])
    xlabel('Time [ms]','fontweight','bold','fontsize',12); 
    ylabel('Angular Velocity [rad/s]','fontweight','bold','fontsize',12)    
end
hold off
legend(BS_Legend, 'Location', 'best')

figure
hold on
for i = 1:4
    % Little Spool
    plot(LS_Time_ms(20:LS_Omega_Max_I(i),i), LS_Omega(20:LS_Omega_Max_I(i),i), 'LineWidth',2,'Color',[colors(i,1),colors(i,2),colors(i,3)])
    xlabel('Time [ms]','fontweight','bold','fontsize',12); 
    ylabel('Angular Velocity [rad/s]','fontweight','bold','fontsize',12)   
end
hold off
legend(LS_Legend, 'Location', 'best')

