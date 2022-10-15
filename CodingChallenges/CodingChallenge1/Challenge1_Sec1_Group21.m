%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE CHALLENGE 1 - 
%
% The purpose of this challenge is to estimate atmospheric pressure in
% Boulder CO using a pressure model and measurements, and compare the two
% through error analysis and statistics. 
%
% To complete the challenge, execute the following steps:
% 1) Load the given dataset
% 2) Extract altitude and pressure data
% 3) Determine standard deviation, variance, mean, and 
%    standard error of the mean of the pressure data
% 4) Using information given about the instrument, find uncertainty associated
%    with altitude measurements
% 5) Use the model to predict pressure measurements at each altitude in the
%    data set, along with propagated uncertainty
% 6) Compare results, discuss, and print answers to the command window. 
% Bonus) Repeat for larger measurement uncertainty in altitude
%
% 1) Avery Gillespie
% 2) Akanksha Nelacanti
% 3) Caleb Bristol
% 4) Ketan Arun Kamat
% 5) Max Bergman
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Housekeeping
clear       % Clear all variables in workspace
close all   % Close all open figure windows
clc         % Clear the command window

%% 1) Load data from given file 

data = readtable("PressureInBoulder.csv");

%% 2) Extract just the altitude and station pressure data columns to meaningfully named variables

AltitudeData = data.Altitude_m_;
PressureData = data.StationPressure_kPa_;

%% 3) Determine Statistics and Error
% the standard deviation, variance, mean, and standard error of the mean (sem) of the pressure data

StdevPressure = std(PressureData);

VarPressure = var(PressureData);

MeanPressure = mean(PressureData);

Sem_Pressure = StdevPressure/sqrt(length(PressureData));

%% 4) Uncertainty
%    The altitude measurements were taken using an instrument that displayed
%    altitude to the nearest tenth of a meter. 

%    What is the associated absolute uncertainty with these measurements? 

AltitudeUncertainty = round(std(AltitudeData)/2,2); % [m]

%% 5) Pressure Predictions 
%    Using the altitude measurements and uncertainty, predict pressure with the follwing model:
%    First, propagate uncertainty BY HAND before calculating uncertainty for each value. 
%    Then check: is it different for each calculation?

%               Model
%    P_est =   P_s * e^(-k*h)               
%    Assume P_s is 101.7 ± 0.4 kPa and k is 1.2*10^(-4) [1/m] 

P_s = 101.7;       % ± 0.4 [kPa]
k = 1.2*10^(-4);   % [1/m]

syms Ps h
P_est = (Ps).*exp(-k*h);
% uncertainty wrt h
dPh = abs(diff(P_est,h)) .* AltitudeUncertainty;
dPh = subs(dPh,h,AltitudeData);
dPh = double(subs(dPh,Ps,P_s));

% uncertainty wrt P_s
dPs = abs(diff(P_est,Ps)) .* 0.4; % P_s uncertainty
dPs = subs(dPs,h,AltitudeData);
dPs = double(subs(dPs,Ps,P_s));

% uncertainty in P_est
P_sig = sqrt(dPh.^2 + dPs.^2);

%% 6) Print Results
%    Display the predicted pressure from the model with it's associated uncertainty and
%    the average pressure with the it's standard error of the mean from the data.

clear P_est h P_s
P_s = 101.7;
P_est = (P_s).*exp(-k.*AltitudeData);

results  = table(P_est,P_sig);
P_data = [num2str(MeanPressure) '  ±  ' num2str(Sem_Pressure) '   kPa'];
disp(results);
disp(P_data);

% Disucss the accuracy of the model and whether or not you think the
% model agrees with the measurements 

disp('Model Discussion: (Write Your Answer Here) ')

%% Bonus
%   Repeat steps 4-6, but assume the altitude measurements were taken on a
%   lower precision instrument that only displayed altitude to nearest 10
%   meters
%   How does this change the results and comparison ? 

altitude_uncertainty_new = round(std(AltitudeData)/2); % [m]

syms Ps h
P_est = (Ps).*exp(-k*h);
% uncertainty wrt h
dPh = abs(diff(P_est,h)) .* AltitudeUncertainty;
dPh = subs(dPh,h,AltitudeData);
dPh = double(subs(dPh,Ps,P_s));

% uncertainty wrt P_s
dPs = abs(diff(P_est,Ps)) .* 0.4; % P_s uncertainty
dPs = subs(dPs,h,AltitudeData);
dPs = double(subs(dPs,Ps,P_s));

% uncertainty in P_est
P_sig = sqrt(dPh.^2 + dPs.^2);

clear P_est
P_est = (P_s).*exp(-k.*AltitudeData);

results  = table(P_est,P_sig);
P_data = [num2str(MeanPressure) '  ±  ' num2str(Sem_Pressure) '   kPa'];
disp(results);
disp(P_data);

