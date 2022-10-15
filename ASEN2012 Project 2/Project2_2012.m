%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 2012 Project 2 - Bottle Rocket  %
%       by Caleb Bristol               %
%       SID: 109599803                 % 
%       Due Date: 12/06/20             %
%       Date Modified: 12/05/20        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Purpose%%
%Plot trajectory and thrust of bottle rocket using ode45
%All variable input will be done through a GUI

%Housekeeping
clear
close all;
clear

%%%%%%%%%%%%%%%%%%%%%%Define Constants%%%%%%%%%%%%%%%%%%%%%%

g = 9.81;                       %[m/s^2]
C_d = 0.8;                      %[Dimensionless] coefficient of discharge
C_D = 0.5;                      %[Dimensionless] coefficient of Drag
rho_air = 0.961;                %[kg/m^3]
rho_h20 = 1000;                 %[kg/m^3]
p_amb = 12.1;                   %[lb/in^2] Ambient Air pressure
p_amb = p_amb * 6894.76;        %[Pa] Convert to consistent units
p_gage_i = 50;                  %[lb/in^2] Initial gage pressure of air in bottle
p_gage_i = p_gage_i * 6894.76;  %[Pa] Better units of pressure
gamma = 1.4;                    %[Dimensionless] Specific Heat ratio for air
R_air = 287;                    %[J/kg*K]
Vol_bot = 0.002;                %[m^3] Volume of Empty Bottle
Vol_h20_i = 0.001;              %[m^3] Iniital volume of water
m_bot_i = 0.15;                 %[kg] Inital mass of empty bottle + fins
T_air_i = 300;                  %[K] Initial temperature of ambient air
D_throat = 2.1;                 %[cm] Diameter of bottle throat
A_throat = pi * 0.25 * (0.01 * D_throat)^2; %[m^2]
D_bot = 10.5;                   %[cm] Diameter of bottle
A_bot = pi * 0.25 * (0.01 * D_bot)^2; %[m^2]
length_stand = 0.5;             %[m] Length of launch stand
theta_0 = 45;                   %[degrees]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%Menu For Altering Variables%%%%%%%%%%%%%%%%%%

options = menu("Run Default Case or Modify Variables","Run Verification Case", "Modify Variables for Test Case");

while options ~= 0
    switch options
        case 1
            options = 0;
        case 2
            modvar = 0;
            while modvar ~= 5
                modvar = menu("Choose Variable to Modify","Air Pressure","Water Volume","Coefficient of Drag","Launch Angle","Exit (Run Test Case)");
                switch modvar
                    case 1
                        fprintf("[Test case pressure is 50 psi] \n");
                        p_gage_i = input("Enter desired initial gage air pressure [psi]: ");
                        p_gage_i = p_gage_i * 6894.76; %Convert to [Pa]
                    case 2
                        fprintf("[Test case Water Volume Fraction is 0.5] \n");
                        frac_h20_i = input("Enter desired initial water volume fraction (from 0 to 1): ");
                        Vol_h20_i = 0.002 * frac_h20_i; %Convert fraction into volume [m^3]
                    case 3
                        fprintf("[Test case Coefficient of Drag is 0.5] \n");
                        C_D = input("Enter desired Coefficient of Drag (from 0.3 to 0.5): ");
                    case 4
                        fprintf("[Test case Launch Angle is 45 degrees] \n");
                        theta_0 = input("Enter desired Launch Angle in degrees (0 to 90): ");
                    case 5
                        break;                      
                end
            end
            break;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%Initial Conditions%%%%%%%%%%%%
v_0 = 0.0;      %[m/s] Initial Velocity
v_x_0 = 0.0;    %with respective x and z components
v_z_0 = 0.0; 
x_0 = 0.0;      %[m]
z_0 = 0.25;     %[m]
theta_0 = theta_0 * (pi / 180); %correct launch angle to [radians]

%Not launchpad length but distance from origin at end of launchpad
launch_stand_x = length_stand * cos(theta_0); %[m] Break into components for r
launch_stand_z = length_stand * sin(theta_0);
launchpad_r = sqrt((launch_stand_x)^2 + (launch_stand_z + z_0)^2); %[m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%Calculate Initial Mass Eq.(12)%%%%%%%%
p_air_i = p_gage_i + p_amb; %[Pa]
Vol_air_i = Vol_bot - Vol_h20_i; %[m^3]
m_rocket_i = m_bot_i + rho_h20 * (Vol_h20_i) + ((p_air_i * Vol_air_i) / (R_air * T_air_i));
m_air_i = (p_air_i * Vol_air_i) / (R_air * T_air_i);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%Set up initial conditions & constants%%%%%%%%%%%%

%constants vector
const = [Vol_air_i p_air_i rho_air C_d C_D A_bot A_throat p_amb rho_h20 g Vol_bot T_air_i R_air m_air_i gamma theta_0 launchpad_r];

%State Vector -> Initial Conditions
X_0 = [x_0 z_0 v_x_0 v_z_0 m_rocket_i m_air_i Vol_air_i];

t = [0 5]; %[s]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%Call ode45%%%%%%%%%%%%%%%%%%%%%%%%%%

[t, X] = ode45(@(t, X)positionfunc(t,X,const), t, X_0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%Modified Function for Creating Thrust Vector%%%%%%%%%

Thrust = zeros(length(t),1);

for i = 1:length(t)
    Thrust(i) = thrustfunc(X(i,:),const);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%Plotting%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Figure 1: Trajectory of rocket
figure(1)
plot(X(:,1),X(:,2)); hold on
xlabel("Distance [m]")
ylabel("Height [m]")
title("Trajectory of Bottle Rocket")
xlim([0 (max(X(:,1)) + 20)])
ylim([0 (max(X(:,2)) + 15)])
grid on
hold off

%Figure 2: Z-velocity of rocket
figure(2)
plot(t,X(:,4)); hold on
xlabel("Time [s]")
ylabel("Z-velocity [m/s]")
title("Bottle Rocket Z-velocity")
xlim([0 5])
grid on
hold off

%Figure 3: Thrust of rocket
figure(3)
plot(t,Thrust); hold on
xlabel("Time [s]")
ylabel("Thrust [N]")
title("Bottle Rocket Thrust vs. Time")
xlim([0 0.5])
ylim([0 inf])
grid on
hold off


%%%%%%%%%%%%Display Maximum x and z values%%%%%%%%%%%%

fprintf("Maximum Height of Rocket [m]: \n")
disp(max(X(:,2)))
fprintf("Maximum Distance of Rocket [m]: \n")
disp(max(X(:,1)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function dXdt = positionfunc(t,X_0,const)
    
%%Read in Constants%%
vol_air_0 = const(1);
p_air_0 = const(2);
rho_air = const(3);
C_d = const(4);
C_D = const(5);
A_bot = const(6);
A_throat = const(7);
p_amb = const(8);
rho_h20 = const(9);
g = const(10);
vol_bot = const(11);
T_air_i = const(12);
R_air = const(13);
m_air_0 = const(14);
gamma = const(15);
launchangle = const(16);
launchpad_r = const(17);

%%Read in State Vector%%
x = X_0(1);
z = X_0(2);
r = sqrt(x^2 + z^2);

v_x = X_0(3);
v_z = X_0(4);

m_rocket = X_0(5);
m_air = X_0(6);
vol_air = X_0(7);

%%Calculate Drag%%
v = sqrt(v_x^2 + v_z^2);
q = 0.5 * rho_air * v^2;
D = q * C_D * A_bot;

if r < launchpad_r  %Define theta while on launch platform (distance traveled < 0.5m)
    theta = launchangle; %[radians]
else
    theta = atan(v_z / v_x);
end

%Calculate p_air for depressurization phase to test if in ballistic phase
p_end = p_air_0 * (vol_air_0 / vol_bot)^gamma;
p_air = p_end * (m_air / m_air_0)^gamma;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%Before Water is Exhausted -> Water Expulsion Phase%%%%%%

if vol_air < vol_bot  
%Phase one: Water Expulsion
    %No air leaving rocket -> no change in air mass
    m_air_dot = 0;
    
    %Define p_air for calculations
    p_air = p_air_0 * (vol_air_0 / vol_air)^gamma;

    %Calculate exit velocity and thus mass flow rate of water
    v_exit = sqrt((2 * (p_air - p_amb))/rho_h20);
    m_h20_dot = C_d * rho_h20 * A_throat * v_exit;
    m_rocket_dot = -m_h20_dot;
    
    %Find thrust vector magnitude
    F = m_h20_dot * v_exit;
    
    %vol_air' = dvdt
    vol_air_dot = C_d * A_throat * v_exit;

    %Calculate Acceleration
    v_x_dot = (F * cos(theta) / m_rocket) - (D * cos(theta) / m_rocket);
    v_z_dot = (F * sin(theta) / m_rocket) - (D * sin(theta) / m_rocket) - g;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%After Water is Exhausted -> Depressurization Phase%%%%%%
    
elseif vol_air >= vol_bot && p_air > p_amb 
%Phase 2: Depressurization
    %Water mass is 0 and no longer changing
    m_h20_dot = 0;
    %Air mass is changing but volume is now simply the bottle volume
    vol_air_dot = 0;
    
    %Define p_end for calculations
    p_end = p_air_0 * (vol_air_0 / vol_bot)^gamma;
    T_end = T_air_i * (vol_air_0 / vol_bot)^(gamma - 1);
    
    %Define pressure, density, and temperature of air at this moment
    p_air = p_end * (m_air / m_air_0)^gamma;
    rho_air = m_air / vol_bot;
    T_air = p_air / (rho_air * R_air);
    
    %Find critical p value
    p_crit = p_air * (2 / (gamma + 1))^(gamma / (gamma - 1));
    

    if p_crit > p_amb
        %Calculate values of p, rho, and T for choked flow
        p_exit = p_crit;
        T_exit = (2 / (gamma + 1)) * T_air;
        rho_exit = p_exit / (R_air * T_exit);
        
        %Calculate V_e for choked flow
        v_exit = sqrt(gamma * R_air * T_exit);
    else
        %Solve for Mach at exit
        M_exit = sqrt(((p_air / p_amb)^((gamma-1)/gamma) - 1) / ((gamma - 1)/2));
        
        %Calculate values of p, rho, and T for unchoked flow
        p_exit = p_amb;
        T_exit = T_air / (1 + ((gamma - 1) / 2) * M_exit^2);
        rho_exit = p_amb / (R_air * T_exit);
        
        %Calculate Exit Velocity for unchoked flow
        v_exit = M_exit * sqrt(gamma * R_air * T_exit);
    end
    
    %Find mass flow rate of air based on (un)choked values found above
    m_air_dot = C_d * rho_exit * A_throat * v_exit;
    m_rocket_dot = -m_air_dot;
    
    %Find thrust vector magnitude
    F = m_air_dot * v_exit + (p_amb - p_exit) * A_throat;
    
    %Convert thrust into acceleration values, add drag and gravity
    v_x_dot = (F * cos(theta) / m_rocket) - (D * cos(theta) / m_rocket);
    v_z_dot = (F * sin(theta) / m_rocket) - (D * sin(theta) / m_rocket) - g;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
%%%%%%%After Bottle is Depressurized -> Ballistic Phase%%%%%%%
    
elseif vol_air >= vol_bot && p_air <= p_amb
%Phase 3: Ballistic Phase
    %No change in air or water -> basically no calculations to be done
    m_h20_dot = 0;
    m_air_dot = 0;
    m_rocket_dot = 0;
    vol_air_dot = 0;

    %Acceleration values without thrust, just drag and gravity
    v_x_dot = - (D * cos(theta) / m_rocket);
    v_z_dot = - (D * sin(theta) / m_rocket) - g;
    
    %Correction for rocket hitting the ground
    if z <= 0
        v_x = 0;
        v_z = 0;
        v_x_dot = 0;
        v_z_dot = 0;
    end
end

%Because mass flows out
m_air_dot = -m_air_dot;
m_h20_dot = -m_h20_dot;

dXdt = [v_x;v_z;v_x_dot;v_z_dot;m_rocket_dot;m_air_dot;vol_air_dot];
end

