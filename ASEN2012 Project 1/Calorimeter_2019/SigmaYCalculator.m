function [PreSigmaY,PostSigmaY] = SigmaYCalculator(Time,Temperature)
[MaxTemp,MaxIndex] = max(Temperature);
MaxTime = Time(MaxIndex);
[MinTemp,MinIndex] = min(Temperature);
MinTime = Time(MinIndex);

%For t = (0,599) there is no object in the calorimeter
%Object inserted at t = 600

%Calorimeter heats up until t = 729 before cooling down again
%goes until t = 1013

%%%%%%%%%%%%%%%Pre Insertion Regression%%%%%%%%%%%%%%%%%%
%Create vectors for time and temp for t = (0,599)
PreTime = Time(1:MinIndex);
PreTemp = Temperature(1:MinIndex);

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
Pre_SigY = sqrt((1/(length(Pre_Y) - length(Pre_x_hat))) * sum(Pre_Deviation .* Pre_Deviation));

PreSigmaY = Pre_SigY;

%%%%%%%%%%%%%%%%%%Post Insertion Regression%%%%%%%%%%%%%%%%%%%
%Create vectors with Post Insertion values t = (729,1013)
PostTime = Time(MaxIndex:714);
PostTemp = Temperature(MaxIndex:714);

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
Post_SigY = sqrt((1/(length(Post_Y) - length(Post_x_hat))) * sum(Post_Deviation .* Post_Deviation));

PostSigmaY = Post_SigY;
end

