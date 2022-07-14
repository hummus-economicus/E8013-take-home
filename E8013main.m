
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

% save surplus matrix
save surplus.mat S_init 


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


%% Question: 8 - simulation 

N_workers = 250;
N_firms = 25;
T = 100; % time periods

% worker matrix entries:
% - worker ID (1st entry)
% - worker productivity (2nd entry)
% - employment status (3rd entry)
% - if employed: Firm ID (4th entry)
% - if employed: Job ID (5th entry)
% - if employed: Firm productivity (6th entry)
% - if employed: Firm job security (high==1, low==2) (7th entry)

M = zeros(N_workers,7,T);
M(:,1,:) = repmat(1:N_workers,1,1,T); % worker IDs
M(:,2,:) = repmat(randi([1 grid_size],1,N_workers),1,1,T); % worker "productivities" / position on grid

% firm matrix entries: 
% - firm ID (1st entry)
% - job ID (2nd entry) 
% - firm productivity  (3rd entry)
% - firm stability: 1=stable  (4th entry)
% - firm match status: 1=matched, 0=unmatched (5th entry)

F = zeros(N_firms,5,T);
F(:,1,:) = repmat(1:N_firms,1,1,T); % firm IDs
F(:,3,:) = repmat(randi([1 grid_size],1,N_firms),1,1,T); % firm "productivities" / position on grid
F(:,4,:) = repmat( rand(N_firms,1)<=phi,1,1,T)*(-1) +2; % firm job security =1 stable, =2 unstable  
F = repmat(kron(F(:,:,1), ones(10,1)) ,1,1,T );
% 
F(:,2,:) = repmat(1:N_firms*10, 1,1,T);% job IDs

S = rand(grid_size,grid_size,2) -0.5; % trial surplus function 
S_plus = max(S,0);     

rng('default')
rand(1)
% 82 = unmatched workers that meet a firm
% ADD SEED 

for t =2:T
    t
   % vector =1 if unmatched worker meets a firm: 
   
   V = (M(:,3,t-1) -1)*(-1).*(rand(N_workers,1)<=lambda);
   nb = sum(V); % sum of unmatched worker that meet a firm 
   
   % pool of unmatched firms 
   sub_firms_unmatched = F(F(:,5,t-1)==0,:,t-1); % select unmatched firms
   % resample randomly firms to simulate random meeting with workers: 
   sub_firms_unmatched_random = sub_firms_unmatched(randperm(size(sub_firms_unmatched, 1)), :);
   % keep only enough firms to match searching workers: 
   sub_firms_meeting = sub_firms_unmatched_random(1:nb,:); %
   
   row = M(V==1,2,t);  % get productivities of workers
   ID = M(V==1,1,t);   % get worker ID
   col = sub_firms_meeting(:,3,1); % get firm productivities
   col2 =  sub_firms_meeting(:,4,1); % get firm job security 
    
   % Among meeting workers, dummy of successful matches: 
   match_success = S(  sub2ind(  size(S), row,col,col2)  )>0; 
   sum(match_success); % 49 sucessful matches 
   % keep IDs of workers with sucessful match:
   ID_matched=ID.*match_success; 
   ID_matched(ID_matched==0) = [];
   % fill out employment status of succesful matches:
   M(ID,3,t) = match_success; 
   % fill out corresponding firm information for each successfully matched
   % worker: 
   M(ID_matched,4:7,t) = sub_firms_meeting(match_success==1,1:4);
   
   % Tag matched firms in firm matrix: 
   F( sub_firms_meeting(match_success==1,2),5,t) = 1;
   
   % Seperate matched workers from last period randomly (depending on job security) 
   if t > 2 
      M(M(:,7,t-1)==1,3,t) =  rand( sum(M(:,7,t-1)==1,1),1)>sigL;
      M(M(:,7,t-1)==2,3,t) =  rand( sum(M(:,7,t-1)==2,1),1)>sigH;
      % If not dissolve, transmit firm info from last period: 
      M(M(:,3,t)==1 & M(:,3,t-1)==1 ,4:7,t) = M(M(:,3,t)==1 & M(:,3,t-1)==1 ,4:7,t-1);
      % transfer status of not dissolved to matched firmds: 
      F( M(  M(:,3,t)==1 , 5,t)  ,5,t) = 1;
   end 
   
end

% Output: 


% define burn-in period:  








