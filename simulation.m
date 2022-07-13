% Question 8

n_f=5000;
n_w=50000;

%
f_mat=zeros(3,n_w);


% start from stationary distribution
% workers: record worker id, x, y (if matched), sigma (if matched), firm id (if matched)
w_mat=zeros(5,n_w);


% simulate 6 years
T=72;
