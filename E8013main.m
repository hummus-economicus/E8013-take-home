
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

% save surplus matrix
save surplus.mat S u_n v_n


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

load surplus.mat
w = equilibrium_wages(S,u_n,v_n,b,alpha,beta,sigL,sigH,phi,lambda,grid_size,tol,tol_out,MaxIt);

% replace negative values by zero

w = max(w,0); 

save surplus.mat S u_n v_n w

% Plot of log-wage for different y and sigma values 
%%
figure; 
subplot(1,2,1)
plot(grid,log(w(:,300,1)))
hold on 
plot(grid,log(w(:,600,1)))
hold on 
plot(grid,log(w(:,900,1)))
title("high job security ")
xlabel('x');
ylabel('log(w(x,y,\sigma_l))');
l=legend('y = 0.2997','y = 0.5994','y = 0.8991');
set(l,'Location','SouthEast');
hold off 
% plot low job security
subplot(1,2,2)
plot(grid,log(w(:,300,2)))
hold on 
plot(grid,log(w(:,600,2)))
hold on 
plot(grid,log(w(:,900,2)))
title("low job security ")
xlabel('x');
ylabel('log(w(x,y,\sigma_h))');
l=legend('y = 0.2997','y = 0.5994','y = 0.8991');
set(l,'Location','SouthEast');
hold off 

%% Question: 8 - simulation 
load surplus.mat

N_workers = 50000;
N_firms = 5000;
T = 12*56; % time periods

[M,F] = simulation(N_workers,N_firms,T,S,u_n,v_n,w,b,alpha,beta,sigL,sigH,phi,lambda,grid_size,tol,tol_out,MaxIt);

%% Question 8b - Regression  

% Keep first month of each year
months = (1:56)*12 - 11;
W = M(:,:,months);
% Eliminate burn-in period of 50 years:
W = W(:,:,51:end);
% reduce to 2 dimensions:
W = reshape(permute(W, [1 3 2]),6*N_workers,7);
% Keep employed workers:  
W(W(:,3)==0,:) = [];

% Create new column for wages 
row = W(:,2);% grid index of worker productivity
col = W(:,6);% grid index of firm productivity 
col2 = W(:,7); % index of job security 
W(:,8) = w(  sub2ind(  size(w), row,col,col2)  ) ; 

% keep variables for regression 
W = W(:, [ 1 2 4 6 7 8]);
W(:,2) = grid(W(:,2));
W(:,4) = grid(W(:,4));
W(:,end) = log(W(:,end));
% 1st column = worker ID
% 2nd column = worker productivity 
% 3rd column = firm ID
% 4th column = firm productivity 
% 5th column = firm job security 
% 6th column = log-wage


T = table(W);
filename = 'regression_table.xlsx';
writetable(T,filename,'Sheet',1,'Range','D1')

%% Question 9: replicate figure VII from Card et al.  

% Eliminate burn-in period of 50 years:
D = M(:,:,50*12+1:end);

% Step 1 
% Identify workers that switched job in 6 year period and held preceeding
% job and new job for two or more years. 

tic;
for i=N_workers:-1:1 % N_workers 
   i
   aux = squeeze(D(i,3,:))';
   
   f = find(diff([0,aux,0]));
   p = f(1:2:end-1);  % Start indices
   y = f(2:2:end)-p;
   
   fi = find(diff([0,y,0]>23));
   pi = fi(1:2:end-1);  % Start indices
   yi = fi(2:2:end)-pi;
   if (isempty(yi) ==0  && yi(1) ~= 2 )
        D(i,:,:) = [];
   elseif isempty(yi)
        D(i,:,:) = []; 
   end
  
end
toc;

save D.mat D
load D

% D = D(1:80,:,:); % TEST <- DELETE LATER
D = D(:,:,[24 49]);

% Create new column for wages 
D = [D , zeros(size(D,1),1,2) ];

row = D(:,2,1);% grid index of worker productivity
col = D(:,6,1);% grid index of firm productivity 
 col2 = D(:,7,1); % index of job security 
D(:,8,1) =   w(  sub2ind(  size(w), row,col,col2 )  )  ; 

row = D(:,2,2);% grid index of worker productivity
col = D(:,6,2);% grid index of firm productivity 
 col2 = D(:,7,2); % index of job security 
D(:,8,2) =   w(  sub2ind(  size(w), row,col,col2 )  )  ; 

% Step 2: Define quartiles based on firm-FE coefficients
% import table as M_firm 

load firmFE.mat 

M_firm=table2array(firmFE);
quants=quantile(M_firm(:,2), [0.25 0.5 0.75]);
M_firm=[M_firm zeros(length(M_firm),4)];
M_firm(:,3)=M_firm(:,2)<=quants(1);
M_firm(:,4)=M_firm(:,2)>quants(1)&M_firm(:,2)<=quants(2);
M_firm(:,5)=M_firm(:,2)>quants(2)&M_firm(:,2)<=quants(3);
M_firm(:,6)=M_firm(:,2)>quants(3);
M_firm(M_firm(:,3)==1,2)=1;
M_firm(M_firm(:,4)==1,2)=2;
M_firm(M_firm(:,5)==1,2)=3;
M_firm(M_firm(:,6)==1,2)=4;
M_firm=M_firm(:,1:2);

% Step 3:  Merge D matrix with firm coefficient using firm ID 

D(D(:,4,1)== 1110 | D(:,4,2)== 1110 ,:,:) = [];
D(D(:,4,1)== 1780 | D(:,4,2)== 1780 ,:,:) = [];
D(D(:,4,1)== 4669 | D(:,4,2)== 4669 ,:,:) = [];

sub_1 = D(:,:,1);
sub_2 =D(:,:,2);
sub_1 = array2table(sub_1);
sub_2 = array2table(sub_2);
M_firm = array2table(M_firm);

M_firm = renamevars(M_firm, ["M_firm1" "M_firm2"], ["firmID", "FE"] );
sub_1 = renamevars(sub_1, ["sub_14" ], ["firmID"] );
sub_2 = renamevars(sub_2, ["sub_24" ], ["firmID"] );

M_pre = join(sub_1,M_firm);
M_post = join(sub_2,M_firm);

% Step 4: keep quartile and wages 
M_pre = M_pre(:,[8 9]);
M_post = M_post(:, [8 9] );
M_pre = table2array(M_pre);
M_post = table2array(M_post);

% 4 - 4
M_pre_4_4 = M_pre(M_pre(:,2)==4 & M_post(:,2)==4 ,:);
M_post_4_4 = M_post(M_pre(:,2)==4 & M_post(:,2)==4,:);
% 4 - 3
M_pre_4_3 = M_pre(M_pre(:,2)==4 & M_post(:,2)==3 ,:);
M_post_4_3 = M_post(M_pre(:,2)==4 & M_post(:,2)==3,:);
% 4 - 2
M_pre_4_2 = M_pre(M_pre(:,2)==4 & M_post(:,2)==2 ,:);
M_post_4_2 = M_post(M_pre(:,2)==4 & M_post(:,2)==2,:);
% 4 - 1 
M_pre_4_1 = M_pre(M_pre(:,2)==4 & M_post(:,2)==1 ,:);
M_post_4_1 = M_post(M_pre(:,2)==4 & M_post(:,2)==1,:);
% 1 - 4
M_pre_1_4 = M_pre(M_pre(:,2)==1 & M_post(:,2)==4 ,:);
M_post_1_4 = M_post(M_pre(:,2)==1 & M_post(:,2)==4,:);
% 1 - 3
M_pre_1_3 = M_pre(M_pre(:,2)==1 & M_post(:,2)==3 ,:);
M_post_1_3 = M_post(M_pre(:,2)==1 & M_post(:,2)==3,:);
% 1 - 2
M_pre_1_2 = M_pre(M_pre(:,2)==1 & M_post(:,2)==2 ,:);
M_post_1_2 = M_post(M_pre(:,2)==1 & M_post(:,2)==2,:);
% 1 - 1
M_pre_1_1 = M_pre(M_pre(:,2)==1 & M_post(:,2)==1 ,:);
M_post_1_1 = M_post(M_pre(:,2)==1 & M_post(:,2)==1,:);

xaxis = -2:1;

figure;
hold all
plot(xaxis,[mean(M_pre_4_4(:,1)) mean(M_pre_4_4(:,1)) mean(M_post_4_4(:,1)) mean(M_post_4_4(:,1))],'--o') 
plot(xaxis,[mean(M_pre_4_3(:,1)) mean(M_pre_4_3(:,1)) mean(M_post_4_3(:,1))  mean(M_post_4_3(:,1))],'--o') 
plot(xaxis,[mean(M_pre_4_2(:,1)) mean(M_pre_4_2(:,1)) mean(M_post_4_2(:,1))  mean(M_post_4_2(:,1))],'--o') 
plot(xaxis,[mean(M_pre_4_1(:,1)) mean(M_pre_4_1(:,1)) mean(M_post_4_1(:,1)) mean(M_post_4_1(:,1))],'--o') 
plot(xaxis,[mean(M_pre_1_4(:,1)) mean(M_pre_1_4(:,1)) mean(M_post_1_4(:,1)) mean(M_post_1_4(:,1))],'--o') 
plot(xaxis,[mean(M_pre_1_3(:,1)) mean(M_pre_1_3(:,1)) mean(M_post_1_3(:,1)) mean(M_post_1_3(:,1))],'--o') 
plot(xaxis,[mean(M_pre_1_2(:,1)) mean(M_pre_1_2(:,1)) mean(M_post_1_2(:,1)) mean(M_post_1_2(:,1))],'--o') 
plot(xaxis,[mean(M_pre_1_1(:,1)) mean(M_pre_1_1(:,1)) mean(M_post_1_1(:,1))  mean(M_post_1_1(:,1))],'--o')
xlabel('Time (0 = first year on new job)');
ylabel('log(w)');
l=legend('4 to 4','4 to 3','4 to 2','4 to 1','1 to 4','1 to 3','1 to 2','1 to 1' );
set(l,'Location','NorthWest');
xticks([-2 -1 0 1])
hold off


%% Question: 12 - redoing for 
sigH = sigH/3;
sigL = sigL/3;
[S,u_n,v_n] = solve_model(b,alpha,beta,sigL,sigH,phi,lambda,grid_size,tol,tol_out,MaxIt);
w = equilibrium_wages(S,u_n,v_n,b,alpha,beta,sigL,sigH,phi,lambda,grid_size,tol,tol_out,MaxIt);
[M,F] = simulation(N_workers,N_firms,T,S,u_n,v_n,w,b,alpha,beta,sigL,sigH,phi,lambda,grid_size,tol,tol_out,MaxIt);
% 12b regression
% Keep first month of each year
months = (1:56)*12 - 11;
W = M(:,:,months);
% Eliminate burn-in period of 50 years:
W = W(:,:,51:end);
% reduce to 2 dimensions:
W = reshape(permute(W, [1 3 2]),6*N_workers,7);
% Keep employed workers:  
W(W(:,3)==0,:) = [];

% Create new column for wages 
row = W(:,2);% grid index of worker productivity
col = W(:,6);% grid index of firm productivity 
col2 = W(:,7); % index of job security 
W(:,8) = w(  sub2ind(  size(w), row,col,col2)  ) ; 

% keep variables for regression 
W = W(:, [ 1 2 4 6 7 8]);
W(:,2) = grid(W(:,2));
W(:,4) = grid(W(:,4));
W(:,end) = log(W(:,end));
% 1st column = worker ID
% 2nd column = worker productivity 
% 3rd column = firm ID
% 4th column = firm productivity 
% 5th column = firm job security 
% 6th column = log-wage


T = table(W);
filename = 'regression_table_newsigma.xlsx';
writetable(T,filename,'Sheet',1,'Range','D1')
