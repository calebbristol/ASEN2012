%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 2012 Project - Calorimeter      %
%       by Caleb Bristol               %
%       SID: 109599803                 % 
%       Date: 10/26/20                 %
%                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Housekeeping
clear
close all;
clear

%%%%%Load Data into Matrices%%%%%
CalorimeterData = load("SampleD");
Time = CalorimeterData(:,1);
TempCal1 = CalorimeterData(:,2);
TempBoilingWater = CalorimeterData(:,3);
TempRoom = CalorimeterData(:,4);
TempCal2 = CalorimeterData(:,5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Here we use the average of the 2 calorimeter measurements for our temp
%We propagate error on account of the averaging later
AvgTempCal = (TempCal1 + TempCal2)./2;

%%%%%Plot Temperature Data%%%%%
figure(1);
plot(Time, AvgTempCal); hold on
title("Calorimeter Temperature Over Time")
xlabel("Time (seconds)")
ylabel("Calorimeter Temperature (C)")
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Here we find the maximum and minimum temperatures and their indices
%The minimum temperature is the approximate time of sample insertion
%The maximum temperature is a point where we assume thermal equilibrium
[MaxTemp,MaxIndex] = max(AvgTempCal);
MaxTime = Time(MaxIndex);
[MinTemp,MinIndex] = min(AvgTempCal);
MinTime = Time(MinIndex);

%For t = (0,599) there is no subject in the calorimeter
%Subject inserted at t = 600

%Calorimeter heats up until t = 729 before cooling down again
%goes until t = 1013

%%%%%%%%%%%%%%%Pre Insertion Regression%%%%%%%%%%%%%%%%%%
%Create vectors for time and temp for t = (0,599)
PreTime = Time(1:MinIndex);
PreTemp = AvgTempCal(1:MinIndex);

% Find number of data points in the vectors
Pre_N = length(PreTime);

% Find linear best fit coefficients A and B
% Create H matrix
Pre_H = [ones(Pre_N,1),PreTime];

% Create y matrix
Pre_Y = PreTemp;

% Create W matrix 
Pre_W = eye(Pre_N);

% Solve for P matrix
Pre_P = (Pre_H' * Pre_W * Pre_H)^-1;

% Solve for x_hat matrix and extract A and B parameters
Pre_x_hat = Pre_P * Pre_H' * Pre_W * Pre_Y;
Pre_A = Pre_x_hat(1);
Pre_B = Pre_x_hat(2);

% extract uncertainty in A and uncertainty in B from P matrix

Pre_Deviation = Pre_Y - (Pre_A + Pre_B .*PreTime);
%Pre_SigY = sqrt((1/(length(Pre_Y) - length(Pre_x_hat))) * sum(Pre_Deviation .* Pre_Deviation));
%This is how we would do it but we propogate error from averaging instead

%Extract Error for temperature set 1 and 2, then average
[PreSigY1,~] = SigmaYCalculator(Time,TempCal1);
[PreSigY2,~] = SigmaYCalculator(Time,TempCal2);
%Here we don't use quadrature because these values are NOT independent
Pre_SigY = PreSigY1 + PreSigY2 * 0.5;

Pre_delta_y (1:length(Pre_Y)) = Pre_SigY;
Pre_Diagonal = 1 ./ (Pre_delta_y .* Pre_delta_y);

%Update weighting matrix
Pre_W = diag(Pre_Diagonal);

%Error Covariance Matrix
Pre_Sigma_xHat = (Pre_H' * Pre_W * Pre_H)^-1;

Pre_A_error = sqrt(Pre_Sigma_xHat(1,1));
Pre_B_error = sqrt(Pre_Sigma_xHat(2,2));

%Predicted Behavior through extrapolation
Pre_Temp_Predicted = Pre_A + Pre_B .* Time;

%%%%%%%%%%%%%%%%%%Post Insertion Regression%%%%%%%%%%%%%%%%%%%
%Create vectors with Post Insertion values t = (729,1013)
PostTime = Time(MaxIndex:714);
PostTemp = AvgTempCal(MaxIndex:714);

% Find number of data points in the vectors
Post_N = length(PostTime);

% Find linear best fit coefficients A and B
% Create H matrix
Post_H = [ones(Post_N,1),PostTime];

% Create y matrix
Post_Y = PostTemp;

% Create W matrix 
Post_W = eye(Post_N);

% Solve for P matrix
Post_P = (Post_H' * Post_W * Post_H)^-1;

% Solve for x_hat matrix and extract A and B parameters
Post_x_hat = Post_P * Post_H' * Post_W * Post_Y;
Post_A = Post_x_hat(1);
Post_B = Post_x_hat(2);

% extract uncertainty in A and uncertainty in B from P matrix

Post_Deviation = Post_Y - (Post_A + Post_B .*PostTime);
%Post_SigY = sqrt((1/(length(Post_Y) - length(Post_x_hat))) * sum(Post_Deviation .* Post_Deviation));
%Same Process as the earlier least-squares fitting
[~,PostSigY1] = SigmaYCalculator(Time,TempCal1);
[~,PostSigY2] = SigmaYCalculator(Time,TempCal2);
Post_SigY = PostSigY1 + PostSigY2 * 0.5;

Post_delta_y (1:length(Post_Y)) = Post_SigY;
Post_Diagonal = 1 ./ (Post_delta_y .* Post_delta_y);

%Create new weighting matrix
Post_W = diag(Post_Diagonal);

%Error Covariance Matrix
Post_Sigma_xHat = (Post_H' * Post_W * Post_H)^-1;

Post_A_error = sqrt(Post_Sigma_xHat(1,1));
Post_B_error = sqrt(Post_Sigma_xHat(2,2));

%Predicted Behavior
Post_Temp_Predicted = Post_A + Post_B .* Time;

%%%%%%%%%%%%%%%%%%%%%%%Plotting%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot of temperature data with inclusion of least-squares fitting lines
figure(2);
plot(Time, AvgTempCal); hold on
title("Calorimeter Temperature Over Time")
xlabel("Time (seconds)")
ylabel("Calorimeter Temperature (C)")
plot(Time,Pre_Temp_Predicted);
plot(Time,Post_Temp_Predicted);
hold off

%%%%%%%%%%%%%%%%%%Final Calculations%%%%%%%%%%%%%%%%%%%%%%
T_avg = (Post_Temp_Predicted(419) + Pre_Temp_Predicted(419)) / 2;
[~,Index_of_Interest] = min(abs(AvgTempCal - T_avg));
%Here we are returned an index of 437, where t = 620
%This is where we will extrapolate to find value of T_2

T_2 = Post_Temp_Predicted(Index_of_Interest);
T_0 = Pre_Temp_Predicted(MinIndex);
%Temperature of sample equivalent to T of boiling water when inserted
T_1 = TempBoilingWater(MinIndex);
Mass_c = 510; %[g]
Mass_s = 88.897; %[g]
C_c = 0.895; %[J/g*C] for AL 6061

%Using Eq. (1)
C_sample = (Mass_c * C_c * (T_2 - T_0)) / (Mass_s * (T_1 - T_2));

%%%%%%%%%%%%%%Final Uncertainty Calculations%%%%%%%%%%%%%%%%

%Sigma for T_2 and T_0 come from the equation
%Sigma_y = sqrt(Sigma_A^2 + x_i^2 * Sigma_B^2)
Sigma_T_2 = sqrt(Post_A_error^2 + Index_of_Interest^2 * Post_B_error^2);
Sigma_T_0 = sqrt(Pre_A_error^2 + MinIndex^2 * Pre_B_error^2);
%Error for T_1 is simply the standard deviation of the boiling water data
%I only included the data up until the point of insertion in this
Sigma_T_1 = std(TempBoilingWater(1:MinIndex));

%SigC^2 = (dC/dT_Max * Sigma_T_Max)^2 + (dC/dT_Min * Sigma_T_Min)^2 + (dC/dT_s * Sigma_T_s)
%See appendix for derivative calculations
dCdT_2 = ((Mass_s * (T_1 - T_2)) * (Mass_c * C_c) - (Mass_c * C_c * (T_2 - T_0)) * -Mass_s) / (Mass_s * (T_1 - T_2))^2;
dCdT_0 = -(Mass_c * C_c) / (Mass_s * (T_1 - T_2));
dCdT_1 = -(Mass_c * C_c) * (T_2 - T_0) / (Mass_s * (T_1 - T_2)^2);

%Final Error of C_Sample
Sigma_C_Sample = sqrt((dCdT_2 * Sigma_T_2)^2 + (dCdT_0 * Sigma_T_0)^2 + (dCdT_1 * Sigma_T_1)^2);

