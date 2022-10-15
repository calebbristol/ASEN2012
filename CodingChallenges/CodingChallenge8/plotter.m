clc 
clear 
close all;

load("result.mat");

x = result.x;
y = result.y;
z = result.z;

figure(1);
plot3(x,y,z,'linewidth',4); hold on
title("Rocket Trajectory")
xlabel("x-axis");
ylabel("y-axis");
zlabel("z-axis");
grid on
hold off