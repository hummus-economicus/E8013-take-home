
clear all; 
clc; 

%% Calibration
b = 0; alpha = 0.5; beta = 0.994; sigL = 0.02; sigH = 0.08;
phi = 0.5; lambda = 0.3; grid_size = 500; 
tol = 10^(-3);


%% Question 4: solve 

% Inner Loop: VFI on Surplus 
grid = linspace(0,1,grid_size);

S_init = zeros(grid_size,grid_size,2); 

% trial with uniform distribution for unemployed and vacancies: 
u = ones(1,grid_size)/grid_size;
v = ones(1,grid_size,2)/grid_size;
v(:,:,1) = v(:,:,1)*phi 
v(:,:,1) = v(:,:,2)*(1-phi) 

S_n = VFI_surplus(S_init,b,alpha,beta,sigL,sigH,phi,lambda,grid_size,u,v,tol);


%% Question 5: Plotting equilibrium matching sets conditional on job security sigma 
up_l = nan(1,grid_size);
low_l = nan(1,grid_size);
up_h = nan(1,grid_size);
low_h = nan(1,grid_size);

% for high job security (sigL)
S_plus = max(S_n,0);
% highest and lowest possible match for given x: 
up_l = max(S_plus(:,:,1),[],2); % max over columns 
low_l= min(S_plus(:,:,1),[],2);
up_h = max(S_plus(:,:,2),[],2); % max over columns 
low_h= min(S_plus(:,:,2),[],2);

figure
subplot(1,2,1)
hold all
plot(grid,low_l)
plot(grid,up_l)
patch([grid fliplr(grid)], [low_l' fliplr(up_l')], 'b')
title("Matching set for high job security")
subplot(1,2,2)
hold all
plot(grid,low_h)
plot(grid,up_h)
patch([grid fliplr(grid)], [low_h' fliplr(up_h')], 'b')
title("Matching set for low job security")
hold off







