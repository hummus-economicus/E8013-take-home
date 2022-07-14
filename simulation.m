function [M,F] = simulation(N_workers,N_firms,T,S,u_n,v_n,w,b,alpha,beta,sigL,sigH,phi,lambda,grid_size,tol,tol_out,MaxIt)

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

end