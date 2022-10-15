%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE CHALLENGE 2 - Monte Carlo Analysis
%
% The purpose of this challenge is to perform a Monte-Carlo analysis on the
% lift generated by an aircraft.  The aircraft has the following characteristics:
%   Wing surface area, S = 80 m^2
%   Lift coefficient, C_L = 0.90 +- 0.03
%
% And is flying under the following conditions
%   Air density, rho = 0.653 kg/m^3
%   Airspeed, V = 100 +- 10 m/s
%
% ---------------------------------------------------------------------------------
%
% To complete the challenge, execute the following steps:
% 1) Sample S, C_L, rho, and V 10,000 times.
% 2) Calculate lift in kilonewtons for each of the 10,000 samplings/simulations.
% 3) Calculate the best estimate and error for lift and report it to the
% command window using appropriate significant figures.
% 4) Plot a histogram of L.
% Bonus 1) Calculate drag in kilonewtons for each of the 10,000
% samplings/simulations.
% Bonus 2) Make a scatterplot of Lift vs Drag.
%
% NOTE: DO NOT change any variable names already present in the code.
% 
% Upload your team's script to Canvas to complete the challenge.
% 
% NAME YOUR FILE AS Challenge2_Sec{section number}_Group{group breakout #}.m 
% ***Section numbers are 1 or 2*** 
% EX File Name: Challenge2_Sec1_Group15.m 
%
%
% 1) Caleb Bristol
% 2) Skylar Clark
% 3) Joshua Pitman
% 4) Madison Paige Davis
% 5) Trevor Reed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Housekeeping 
% (Please don't "clear all" or "clearvars", it makes grading difficult)
close all   % Close all open figure windows
clc         % Clear the command window

%% 1) Sample S, C_L, rho, and V 10,000 times
% (i.e. the S variable should contain 10000 samples of the wing surface area)
N = 1e04;
S = 80 .* ones(N,1); %m^2   %despite our assumption of exact values, S and rho given sample vectors N long
C_L = 0.9 + (0.06 .* randn(N,1) - 0.03);    %error simulated by random normal distribution between margin values
rho = 0.653 .* ones(N,1);    %kg/m^3
V = 100 + (20 .* randn(N,1) - 10);

%% 2) Calculate lift in kilonewtons for each of the 10,000 samplings/simulations.
% Given that the equation for lift is:
    %L = 0.5 .* rho .* V.^2 .* C_L .* S; (Newtons)

L = 0.5 .* rho .* V.^2 .* C_L .* S; %(Newtons)

%% 3) Calculate the best estimate and error for lift
% Report it to the command window using appropriate significant figures.
L_best = mean(L);
L_best = round(L_best,3,'significant');     %calculate avg value for best, and round to 3 sig figs
    fprintf("The best value for Lift is %d \n", L_best);
L_err  = std(L);
L_err = round(L_err,2,'significant');   %rounded to 2 sig figs to be same order of magnitude as L_best
    fprintf("The standard deviation of Lift is %d \n", L_err);

%% 4) Plot a histogram (use the "histogram" command) of L with 30 bins.  
% Add annotations and labels for style points!
figure(1);
histogram(L,30,'FaceColor','r'); hold on
    xlabel("Lift (N)");
    ylabel("#occurences");
    title("Monte Carlo Histogram - Lift");
    ann = ["Best value for Lift (N): ", num2str(round(L_best)), "Standard deviation of Lift (N): ", num2str(round(L_err))];
    annotation('textbox',[.6 .5 .3 .3], 'String', ann ,'FitBoxToText','on')
    hold off;

%% Bonus 1) Calculate drag in kilonewtons 
% For each of the 10,000 samplings/simulations, given that the equation for drag is:
%       D = 0.5 * rho * V^2 * C_D * S (Newtons)
% and that C_D = 0.070 +- 0.005
C_D = 0.070 + (0.010 .* randn(N,1) - 0.005);
D = 0.5 .* rho .* V.^2 .* C_D .* S; %(Newtons)

%% Bonus 2) Make a scatterplot of Lift vs Drag.
figure(2);  
scatter(L, D, 5, 'magenta', 'filled');   hold on
    xlabel("Lift (N)");
    ylabel("Drag (N)");
    title("Scatterplot: Lift vs Drag");
hold off;

figure(3); hold on    %bonus linear regression model because it seemed fitting
    lobf = fitlm(L,D);  
    plot(lobf);
    xlabel("Lift (N)");
    ylabel("Drag (N)");
    title("Linear Regression Model: Lift vs Drag");
    hold off;
%
% Think about the following (no work to do):
%     - Why do you think the points are spread into an ellipse and not a
%     circle?
%     - What is the significance of the general trend/slope of the data?
%     - How could this sort of analysis be useful when dealing with more
%     complicated systems and equations?





















