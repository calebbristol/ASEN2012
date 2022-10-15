function T = thrustfunc(X_0,const)
    
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
    F = 0;
    
    %Correction for rocket hitting the ground
    if z <= 0
        v_x = 0;
        v_y = 0;
        v_x_dot = 0;
        v_z_dot = 0;
    end
end

%Because mass flows out
m_air_dot = -m_air_dot;
m_h20_dot = -m_h20_dot;

T = F;
end

