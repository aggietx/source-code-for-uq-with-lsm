% clear;
close all;
% rng(10);
Linewid=2;fonsize=18;
% Generating a true signal of the SPEKF-M model
T = 500;
dt = 0.01;
N = round(T/dt);
x_truth = zeros(1,N);
y_truth = zeros(1,N);
z_truth = zeros(1,N);
x_truth(1,1)=truex;
z_truth(1,1)=truez;
y_truth(1,1)=truey;
Tacf=15;range=floor(50/dt):floor(150/dt);
rho = 28;
sigma = 10;
beta = 8/3;

% sigma_x = 5;
% sigma_y = 5;
% sigma_z = 5;


% sigma_x = 10;
% sigma_y = 1;
% sigma_z = 1;

sigma_x = 0;
sigma_y = 0;
sigma_z = 0;

for i = 2:N
    x_truth(:,i) = x_truth(:,i-1) + sigma * (y_truth(:,i-1) - x_truth(:,i-1)) * dt + sigma_x * sqrt(dt) * randn; 
    y_truth(:,i) = y_truth(:,i-1) + (x_truth(:,i-1) * (rho - z_truth(:,i-1)) - y_truth(:,i-1)) * dt + sigma_y * sqrt(dt) * randn;
    z_truth(:,i) = z_truth(:,i-1) + (x_truth(:,i-1) * y_truth(:,i-1) - beta * z_truth(:,i-1)) * dt + sigma_z * sqrt(dt) * randn;
end