%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE CHALLENGE 5 - Template Script
%
% The purpose of this challenge is to predict whether or not the Boulder
% Reservior will have to close due to a major leak.
%
% To complete the challenge, execute the following steps:
% Part 1:
% 1) Read in the data file
% 2) Set values to any constants
% 3) Perform a trapazoid integration on the data w/r.t. x
% 4) Perform a simpson's 1/3 integration on the data w/r.t. x
% 5) Display which volume measurement is more accurate and why
%
% Part 2:
% 1) Define which delta t will be used in the Euler integration
% 2) Set values to any constants and initial conditions
% 3) Propagate h with t using Euler integration
% 4) Repeat steps 1-4 with different delta t values
% 5) Display which delta t gives a more accurate result and why.
% 
%
% NOTE: DO NOT change any variable names already present in the code.
% 
% Upload your team's script to Gradescope to complete the challenge.
% 
% NAME YOUR FILE AS Challenge5_Sec1_Group21.m 
% ***Section numbers are 1 or 2*** 
% EX File Name: Challenge5_Sec1_Group15.m 
%
%
% 1) Andrew Nichols
% 2) Sean Svihla
% 3) Caleb Bristol
% 4) Ruize Liu
% 5) Chris Lolkema
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Housekeeping
% don't "clear variables", it makes things easier to grade
close all;   % Close all open figure windows
clc;         % Clear the command window
%% Part 1

%% Set up
data = readtable("depth_data.csv"); % read in .csv
x = data.x; % [ft]
d = data.d; % [ft]
L = 4836; % length of reservior [ft]
%% Trapazoid - Calculate Volume
n = length(data.d)-1;
delta_x = (data.x(end)-data.x(1)) / 30;
Area_trap = 0;
for i = 1:n
    Area_trap = Area_trap + (delta_x/2 * (data.d(i) + data.d(i+1))); 
end
Vol_trap = L*Area_trap;% [ft^3]


%% Simpson 1/3 - Calculate Volume
Area_simp = 0;
for i = 2:2:n
    Area_simp = Area_simp + delta_x/3 * (4*data.d(i) + data.d(i+1) + data.d(i-1)); 
end
Vol_simp = L*Area_simp;% [ft^3]


%% Part 2

%% Set up
del_t = [.5 1 4 7]; % various delta t values to test [days]
figure(1) % create figure
h0 = 20; % initial depth
alpha = 1.5e6; % relating volume out per day to depth [ft^2/day]
dV_in = 2e7; % volume in rate per day
for j = 1:4
    t = []; % allocate time vector [days]
    t = 0:del_t(j):200; % allocate time vector [days]
    h = []; % allocate depth vector [ft]
    h(1)= h0; % set initial value in h vector
    for i = 1:(length(t)-1) % Euler method
        dhdt = get_dhdt(h(i),L,alpha,dV_in); % get dh/dt at this depth
        h(i+1) = h(i) + dhdt * del_t(j); %compute next depth value
    end
% plot results
fprintf("The water level evens out at a depth of: %.3f \n",h(end));
plot(t,h);
hold on
end
hold off
% labels for plot
legend(".5 days", "1 days", "4 Days", "7 days")
title("Height vs Time for Boulder Reservoir")

fprintf(".5 Days is the most accurate, however it is a waste of processing power unless you have a NASA super computer")
fprintf("\n They all eventually even out at the same depth, the difference between the different del_t is just what day the reservoir evens out.\n")





