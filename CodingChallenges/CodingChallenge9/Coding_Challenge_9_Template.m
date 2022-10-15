%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE CHALLENGE 9 - Guided Template Script
%
% The purpose of this challenge is to propagate an orbit in a two body
% system for one period, and to plot it's specific energy over time.
%
% To complete the challenge, execute the following steps:
% 1) Set an initial condition vector
% 2) Propagate for exactly period of the orbit
% 3) Calculate the specific energy of the s/c vs. time
% 4) Plot the trajectory, include points for where the trajectory starts,
% ends, and the where the Earth is.
% 5) Plot the change in specific energy vs. time
%
% NOTE: DO NOT change any variable names already present in the code.
% 
% Upload your team's script to Gradescope when complete.
% 
% NAME YOUR FILE AS Challenge9_Sec{section number}_Group{group breakout #}.m 
% ***Section numbers are 1 or 2*** 
% EX File Name: Challenge9_Sec1_Group15.m 
%
% STUDENT TEAMMATES
% 1) Caleb Bristol
% 2) Andrew Nichols
% 3) Atkin Arnstein
% 4) David Hightower
% 5) Mrityunjaya Lala
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Housekeeping
clear variables; close all; clc;


%% Set up
mu = 398600.4415; % GM of Earth [km^3/s^2]
r = [10000; 0; 1000]; % initial r vetor [km]
v = [0; 7.5574; 0]; % initial v vetor [km/s]
a = -mu * (abs(v(1)^2 + v(2)^2 + v(3)^2) - 2 * (mu / sqrt(r(1)^2 + r(2)^2 + r(3)^2)))^-1; % calculating a [km]
T = 2 * pi * sqrt(a^3 / mu); % calculating T [s]
IC = [r; v]; % initial condition vector
t = [0 T]; % time domain [s]

%% Propagate w/ ode45
[t, orbit] = ode45(@(t,orbit)fungi(t,orbit,mu), t, IC);

%% Calculate specific energy
energy = abs(orbit(:,4).^2 + orbit(:,5).^2 + orbit(:,6).^2) - 2 .* (mu ./ sqrt(orbit(:,1).^2 + orbit(:,2).^2 + orbit(:,3).^2));

%% Plotting
figure(1)
plot3(r(1), r(2), r(3),'*'); % plot starting point
hold on; grid minor;
plot3(0, 0, 0,'*'); % plot earth
xlabel("x-axis [km]"); 
ylabel("y-axis [km]");
zlabel("z-axis [km]");
title("Orbit of Satellite about Earth");
plot3(orbit(:,1), orbit(:,2), orbit(:,3)); % plot trajectory
plot3(orbit(end,1), orbit(end,2), orbit(end,3),'*'); % plot ending point
legend("Starting Point", "Earth", "Orbit", "Ending Point");


figure(2)
plot(t, energy); %plot specific energy vs. time
grid minor;
xlabel("Time (s)");
ylabel("Specific Energy [J/kg]");
title("Specific Energy of Orbiting Body");
 
function drdt = fungi(t,IC,mu)
r = sqrt(IC(1)^2 + IC(2)^2 + IC(3)^2);
%vx' = -GM/r^3 * rx
%vy' = -GM/r^3 * ry
%vz' = -Gm/r^3 * rz
vx = IC(4);
vy = IC(5);
vz = IC(6);

vx_dot = -mu / r^3 * IC(1);
vy_dot = -mu / r^3 * IC(2);
vz_dot = -mu / r^3 * IC(3);

drdt = [vx;vy;vz;vx_dot;vy_dot;vz_dot];
end