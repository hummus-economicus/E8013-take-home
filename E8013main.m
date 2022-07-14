
clear all; 
clc; 

%% Parameters
b = 0; alpha = 0.5; beta = 0.994; sigL = 0.02; sigH = 0.08;
phi = 0.5; lambda = 0.3; grid_size = 1000; 
tol = 10^(-3);
tol_out = 10^(-5); % tolerance value for outer loop
MaxIt = 10^4; 


%% Question 4: solve 

grid = linspace(grid_size^(-1),1-grid_size^(-1),grid_size);
[S,u_n,v_n] = solve_model(b,alpha,beta,sigL,sigH,phi,lambda,grid_size,tol,tol_out,MaxIt);


%% Question 5: Plotting equilibrium matching sets conditional on job security sigma 
up_l = nan(1,grid_size);
low_l = nan(1,grid_size);
up_h = nan(1,grid_size);
low_h = nan(1,grid_size);

S_plus = max(S,0);

for x = 1:grid_size
    if isempty(find(S_plus(x,:,1),1,'last')); up_l(x) = NaN; else; up_l(x) = grid(find(S_plus(x,:,1),1,'last')); end
    if isempty(find(S_plus(x,:,1),1)); low_l(x) = NaN; else; low_l(x) = grid(find(S_plus(x,:,1),1)); end

    if isempty(find(S_plus(x,:,2),1,'last')); up_h(x) = NaN; else; up_h(x) = grid(find(S_plus(x,:,2),1,'last')); end
    if isempty(find(S_plus(x,:,2),1)); low_h(x) = NaN; else; low_h(x) = grid(find(S_plus(x,:,2),1)); end
end 
grid_l = grid(~isnan(low_l));
grid_h = grid(~isnan(low_h));

figure
subplot(1,2,1)
hold all
plot(grid,low_l,'Color',[0 0 0])
plot(grid,up_l,'Color',[0 0 0])
patch([grid_l fliplr(grid_l)]', [up_l(~isnan(low_l)) fliplr(low_l(~isnan(low_l)))]', [.5 .5 .5])
title("high job security")
subplot(1,2,2)
hold all
plot(grid,low_h,'Color',[0 0 0])
plot(grid,up_h,'Color',[0 0 0])
patch([grid_h fliplr(grid_h)]', [up_h(~isnan(low_h)) fliplr(low_h(~isnan(low_h)))]', [.5 .5 .5])
title("low job security")
hold off


%% Question 7: equilibrium wages

w = equilibrium_wages(S,u_n,v_n,b,alpha,beta,sigL,sigH,phi,lambda,grid_size,tol,tol_out,MaxIt);

% Plot of log-wage for different y and sigma values 

figure; 
subplot(1,2,1)
plot(grid(S(:,1,1)>=0),log(w(S(:,1,1)>=0,1,1)))
hold on 
plot(grid(S(:,200,1)>=0),log(w(S(:,200,1)>=0,200,1)))
hold on 
plot(grid(S(:,400,1)>=0),log(w(S(:,400,1)>=0,400,1)))
title("high job security ")
xlabel('x');
ylabel('\log(w(x,y,\sigma_l))');
l=legend('y = 0','y = 0.3988','y = 0.7996');
set(l,'Location','SouthEast');
hold off 
subplot(1,2,2)
plot(grid(S(:,1,2)))>=0),log(w(S(:,1,2)))>=0,1,2)))
hold on 
plot(grid(S(:,200,2)>=0),log(w(S(:,200,2)>=0,200,2)))
hold on 
plot(grid(S(:,400,2)>=0),log(w(S(:,400,2)>=0,400,2)))
title("low job security ")
xlabel('x');
ylabel('\log(w(x,y,\sigma_h))');
l=legend('y = 0','y = 0.3988','y = 0.7996');
set(l,'Location','SouthEast');
hold off 


%% Question: 8 - simulation 

N_workers = 250;
N_firms = 25;
T = 100; % time periods

[M,F] = simulation(N_workers,N_firms,T,S,u_n,v_n,w,b,alpha,beta,sigL,sigH,phi,lambda,grid_size,tol,tol_out,MaxIt);




%% Question: 12 - redoing for 
sigH = sigH/3;
sigL = sigL/3;
[S,u_n,v_n] = solve_model(b,alpha,beta,sigL,sigH,phi,lambda,grid_size,tol,tol_out,MaxIt);
w = equilibrium_wages(S,u_n,v_n,b,alpha,beta,sigL,sigH,phi,lambda,grid_size,tol,tol_out,MaxIt);
[M,F] = simulation(N_workers,N_firms,T,S,u_n,v_n,w,b,alpha,beta,sigL,sigH,phi,lambda,grid_size,tol,tol_out,MaxIt);