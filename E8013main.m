
clear all; 
clc; 

%% Calibration
b = 0; alpha = 0.5; beta = 0.994; sigL = 0.02; sigH = 0.08;
phi = 0.5; lambda = 0.3; grid_size = 500; 
tol = 10^(-3);

%% Inner Loop: VFI on Surplus 

% set intial surplus to production value f(x,y) = xy
S_init = repmat(grid,grid_size,1).*repmat(grid',1,grid_size);

% trial with uniform distribution for unemployed and vacancies: 
u = ones(1,grid_size)/grid_size;
v = ones(1,grid_size,2)/grid_size;
S_n = VFI_surplus(S_init,b,alpha,beta,sigL,sigH,phi,lambda,grid_size,u,v,tol);



%% Plotting equilibrium matching sets conditional on job security sigma 










