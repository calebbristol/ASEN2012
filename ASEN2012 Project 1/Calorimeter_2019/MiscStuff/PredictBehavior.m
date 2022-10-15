function [Pre_Temp_Predicted,Pre_A,Pre_A_error,Pre_B,Pre_B_error,Post_Temp_Predicted,Post_A,Post_A_error,Post_B,Post_B_error] = PredictBehavior(Time,temperature)

[MaxTemp,MaxIndex] = max(temperature);
MaxTime = Time(MaxIndex);
[MinTemp,MinIndex] = min(temperature);
MinTime = Time(MinIndex);

%For t = (0,599) there is no object in the calorimeter
%Object inserted at t = 600

%Calorimeter heats up until t = 729 before cooling down again
%goes until t = 1013

%%%%%%%%%%%%%%%Pre Insertion Regression%%%%%%%%%%%%%%%%%%
%Create vectors for time and temp for t = (0,599)
PreTime = Time(1:MinIndex);
PreTemp = temperature(1:MinIndex);

% Find number of data points in the vectors
Pre_N = length(PreTime);

% Find linear best fit coefficients A and B
% Create H matrix
Pre_H = [ones(Pre_N,1),PreTime];

% Create y matrix
Pre_Y = PreTemp;

% Create W matrix (hint: type <help diag> in command line)
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

Pre_delta_y (1:length(Pre_Y)) = Pre_SigY;
Pre_Diagonal = 1 ./ (Pre_delta_y .* Pre_delta_y);

Pre_W = diag(Pre_Diagonal);

%Error Covariance Matrix
Pre_Sigma_xHat = (Pre_H' * Pre_W * Pre_H)^-1;

Pre_A_error = sqrt(Pre_Sigma_xHat(1,1));
Pre_B_error = sqrt(Pre_Sigma_xHat(2,2));

%Predicted Behavior
Pre_Temp_Predicted = Pre_A + Pre_B .* Time;

%%%%%%%%%%%%%%%%%%Post Insertion Regression%%%%%%%%%%%%%%%%%%%
%Create vectors with Post Insertion values t = (729,1013)
PostTime = Time(MaxIndex:714);
PostTemp = temperature(MaxIndex:714);

% Find number of data points in the vectors
Post_N = length(PostTime);

% Find linear best fit coefficients A and B
% Create H matrix
Post_H = [ones(Post_N,1),PostTime];

% Create y matrix
Post_Y = PostTemp;

% Create W matrix (hint: type <help diag> in command line)
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

Post_delta_y (1:length(Post_Y)) = Post_SigY;
Post_Diagonal = 1 ./ (Post_delta_y .* Post_delta_y);

Post_W = diag(Post_Diagonal);

%Error Covariance Matrix
Post_Sigma_xHat = (Post_H' * Post_W * Post_H)^-1;

Post_A_error = sqrt(Post_Sigma_xHat(1,1));
Post_B_error = sqrt(Post_Sigma_xHat(2,2));

%Predicted Behavior
Post_Temp_Predicted = Post_A + Post_B .* Time;
end

