
clear all; 
clc; 

%% Calibration
b = 0; alpha = 0.5; beta = 0.994; sigL = 0.02; sigH = 0.08;
phi = 0.5; lambda = 0.3; grid_size = 500; 
tol = 10^(-3);
tol_out = 10^(-15); % tolerance value for outer loop
MaxIt = 10^4; 


%% Question 4: solve 

grid = linspace(0,1,grid_size);

% Initial guess: 
S_init = zeros(grid_size,grid_size,2); 
u_init = ones(1,grid_size)/grid_size;
v_init = ones(1,grid_size,2)/grid_size;
v_init(:,:,1) = v_init(:,:,1)*phi ;
v_init(:,:,1) = v_init(:,:,2)*(1-phi) ;

% Inner Loop: VFI on Surplus 
u_n = u_init; 
v_n = v_init;

v_n1 = nan(grid_size,1,2);
u_n1 = nan(grid_size,1);

check=1;
it = 0
while (check > tol_out & it < MaxIt)
    it = it + 1; 
    S_n1 = VFI_surplus(S_init,b,alpha,beta,sigL,sigH,phi,lambda,grid_size,u_n,v_n,tol,MaxIt);

% update distributions 
    
    for y = 1:grid_size 
       auxL = S_n1(:,y,1)>=0./( S_n1(:,y,1)<0 +sigL*(S_n1(:,y,1)>=0) ) 
       auxH = S_n1(:,y,2)>=0./( S_n1(:,y,2)<0 +sigH*(S_n1(:,y,2)>=0) ) 

       v_n1(:,:,1) = phi/(1 + u_n*lambda*auxL )
       v_n1(:,:,2) = (1-phi)/(1 + u_n*lambda*auxH )
    end 
    for x = 1:grid_size
        auxL = S_n1(x,:,1)>=0./( S_n1(x,:,1)<0 +sigL*(S_n1(x,:,1)>=0) ) 
        auxH = S_n1(x,:,2)>=0./( S_n1(x,:,2)<0 +sigH*(S_n1(x,:,2)>=0) ) 
        u_n1(x,:) = 1/(1 + phi*v_n(:,:,1)*lambda*auxL' + (1-phi)*v_n(:,:,2)*lambda*auxH' )
    end
    
    check_v = max(  abs( v_n1-v_n ) , [], 1 );
    check_u = max(  abs( u_n1-u_n ) , [], 1 );

    check = max([check_v,check_u]);

    S_init = S_n1; 
    u_n = u_n1;
    v_n = v_n1;
    
end 


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



