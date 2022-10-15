%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE CHALLENGE 7 - Template Script
%
% The purpose of this challenge is to estimate the velocity and kenitic
% energy profile of a falling object. 
%
% To complete the challenge, execute the following steps:
% 1) Set an initial condition velocity
% 2) Set values for constants
% 3) Propagate freefall w/ drag for 20 seconds
% 4) Plot the velocity vs. time
% 5) Calculate the change kinetic energy vs. time
% 6) Plot the change in kinetic energy vs. time
%
% NOTE: DO NOT change any variable names already present in the code.
% 
% Upload your team's script to Gradescope when complete.
% 
% NAME YOUR FILE AS Challenge7_Sec{section number}_Group{group breakout #}.m 
% ***Section numbers are 1 or 2*** 
% EX File Name: Challenge7_Sec1_Group15.m 
%
% STUDENT TEAMMATES
% 1) Sean Mccluskey
% 2) Caleb Bristol
% 3) Hattie Rice
% 4) Tyler Schwinck
% 5) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Housekeeping
clear variables; close all; clc;

%% Set up
m = 0.3; % [kg]
g = 9.81; % [m/s^2]
rho = 1.225; % [kg/m^3]
Cd = 1.2; % coefficient of drag
A = 0.0046; % [m^2]
v0 = 0; % [m/s]

%% Propagate with ode45
t0 = 0;
tf = 20;
dt = 1;

[t,V] = ode45(@(t,V)forceFunc(t,V,Cd,rho,A,m,g), [t0:dt:tf], v0);


%% Plot Velocity vs. Time
figure(1)
plot(t,V); hold on
xlabel("Time [sec]");
ylabel("Velocity [m/s]");
title("Velocity vs. Time");
hold off

%% Calculate Kinetic Energy 
KE = 0.5 .* m .* V.^2;

%% Plot Kinetic Energy vs. Time
figure(2)
plot(t,KE); hold on
xlabel("Time [sec]");
ylabel("Kinetic Energy [J]");
title("Kinetic Energy vs. Time");
hold off
