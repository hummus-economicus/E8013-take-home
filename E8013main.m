
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

% Initial guess: 
S_init = zeros(grid_size,grid_size,2); 
u_init = ones(1,grid_size)/grid_size;
v_init = ones(1,grid_size,2)/grid_size;
v_init(:,:,1) = v_init(:,:,1)*phi ;
v_init(:,:,1) = v_init(:,:,2)*(1-phi) ;

% Inner Loop: VFI on Surplus 
u_n = u_init; 
v_n = v_init;

u_n1 = nan(1,grid_size);
v_n1 = nan(1,grid_size,2);

check=1;
it = 0
while (check > tol_out)
    if it == MaxIt
        disp('Failed to converge!')
        break
    end 

    S_n1 = VFI_surplus(S_init,b,alpha,beta,sigL,sigH,phi,lambda,grid_size,u_n,v_n,tol,MaxIt);

% update distributions 
    
    v_n1(:,:,1) = (( 1 + lambda*u_n*( (S_n1(:,:,1)>=0)./( (S_n1(:,:,1)<0) + sigL*(S_n1(:,:,1)>=0) ) ) )/phi).^(-1);
    v_n1(:,:,2) = (( 1 + lambda*u_n*( (S_n1(:,:,2)>=0)./( (S_n1(:,:,2)<0) + sigH*(S_n1(:,:,2)>=0) ) ) )/(1-phi)).^(-1);

    auxL = (S_n1(:,:,1)>=0)./( (S_n1(:,:,1)<0) + sigL*(S_n1(:,:,1)>=0) );
    auxH = (S_n1(:,:,2)>=0)./( (S_n1(:,:,2)<0) + sigH*(S_n1(:,:,2)>=0) );
    u_n1 = ( 1 + phi*lambda*v_n(:,:,1)*auxL' + (1-phi)*lambda*v_n(:,:,2)*auxH' ).^(-1);
    
    check_v = max(  abs( v_n1-v_n ) , [], 'all' );
    check_u = max(  abs( u_n1-u_n ) , [], 'all' );

    check = max([check_v,check_u])

    S_init = S_n1; 
    u_n = u_n1;
    v_n = v_n1;

    it = it + 1
    
end; if check < tol_out; disp('converged.'); end 


%% Question 5: Plotting equilibrium matching sets conditional on job security sigma 
up_l = nan(1,grid_size);
low_l = nan(1,grid_size);
up_h = nan(1,grid_size);
low_h = nan(1,grid_size);

S_plus = max(S_init,0);

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
title("Matching set for high job security")
subplot(1,2,2)
hold all
plot(grid,low_h,'Color',[0 0 0])
plot(grid,up_h,'Color',[0 0 0])
patch([grid_h fliplr(grid_h)]', [up_h(~isnan(low_h)) fliplr(low_h(~isnan(low_h)))]', [.5 .5 .5])
title("Matching set for low job security")
hold off


%% Question 7: equilibrium wages

w = nan(grid_size,grid_size,2);
V = phi*sum(v_n(:,:,1),'all') + (1-phi)*sum(v_n(:,:,2),'all');

for x = 1:grid_size
   for y = 1:grid_size
       w(x,y,1) = (1- beta*(1-sigL))*alpha*Sn1(x,y,1) + b ...
                   + alpha*beta*lambda*( phi*S_plus(x,:,1)*v_n(:,:,1)'/V ...
                        + (1-phi)*S_plus(x,:,2)*v_n(:,:,2)'/V) 
       w(x,y,2) = (1- beta*(1-sigH))*alpha*Sn1(x,y,2) + b ...
                   + alpha*beta*lambda*( phi*S_plus(x,:,1)*v_n(:,:,1)'/V ...
                        + (1-phi)*S_plus(x,:,2)*v_n(:,:,2)'/V) 
   end
end 

% Plot of log-wage for different y and sigma values 

figure; 
subplot(1,2,1)
plot(grid,log(w(x,1,1)))
hold on 
plot(grid,log(w(x,200,1)))
hold on 
plot(grid,log(w(x,400,1)))
title("Equilibrium log-wage for high job security ")
xlabel('x');
ylabel('log(w(x,y,sigma_l))');
l=legend('y = 0','y = 0.3988','y = 0.7996');
set(l,'Location','SouthEast');
hold off 
subplot(1,2,2)
plot(grid,log(w(x,1,2)))
hold on 
plot(grid,log(w(x,200,2)))
hold on 
plot(grid,log(w(x,400,2)))
title("Equilibrium log-wage for low job security ")
xlabel('x');
ylabel('log(w(x,y,sigma_h))');
l=legend('y = 0','y = 0.3988','y = 0.7996');
set(l,'Location','SouthEast');
hold off 

