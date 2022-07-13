function S_n = VFI_surplus(S_init,b,alpha,beta,sigL,sigH,phi,lambda,grid_size,u,v,tol,MaxIt)
% -------------------------------------------------------------------------
% This function performs a value function ...
% 
% Inputs:
% - S_init: initial guess for surplus function 
% - ... other parameters from exercise 
% - u: density of unemployed (1xgrid_size) vector
% - v: density of vacancies  (1xgrid_sizex2) vector
% - tol: tolerance value for convergence 
% - MaxIt: Maximum number of iterations 
% -------------------------------------------------------------------------

% Creating equidistant grid for productivities x,y
grid = linspace(0,1,grid_size);

flow = nan(grid_size, grid_size);
S_plus = nan(grid_size, grid_size,2);
S_n1 = nan(grid_size, grid_size, 2);

% x in rows
% y in columns 
flow = repmat(grid,grid_size,1).*repmat(grid',1,grid_size)-b;

S_n = S_init;

check_S=1;
it = 0
while (check_S > tol & it < MaxIt)
    it = it + 1; 
    % define S+ = max(S(x,y,sig),0):
    S_plus = max(S_n,0);
    % define total number of unemployed and vacancies:
    U = sum(u); 
    V = phi*sum(v(:,:,1),'all') + (1-phi)*sum(v(:,:,2),'all');
    
for x = 1:grid_size
    for y = 1:grid_size 
        S_n1(x,y,1) = flow(x,y)  ...
                        - (1-alpha)*beta*lambda*S_plus(:,y,1)'*u'/U ...
                        - alpha*beta*lambda*( phi*S_plus(x,:,1)*v(:,:,1)'/V ...
                        + (1-phi)*S_plus(x,:,2)*v(:,:,1)'/V) ;
        S_n1(x,y,2) = flow(x,y)  ...
                        - (1-alpha)*beta*lambda*S_plus(:,y,1)'*u'/U ...
                        - alpha*beta*lambda*( phi*S_plus(x,:,1)*v(:,:,1)'/V ...
                        + (1-phi)*S_plus(x,:,2)*v(:,:,1)'/V) ;
        S_n1(x,y,1) = (1-(1-sigL)*beta)^(-1)*S_n1(x,y,1) ;
        S_n1(x,y,2) = (1-(1-sigH)*beta)^(-1)*S_n1(x,y,2) ;
    end 
end
    
    check_S_h = max(  abs( S_n1(:,:,1)-S_n(:,:,1) ) , [], 1 );
    check_S_l = max(  abs( S_n1(:,:,2)-S_n(:,:,2) ) , [], 1 );
    check_S = max(check_S_h,check_S_l);
    
    S_n = S_n1;
    
end

end
