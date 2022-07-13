function S_n = VFI_surplus(S_init,b,alpha,beta,sigL,sigH,phi,lambda,grid_size,u_n,v_n,tol,MaxIt)
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

S_plus = nan(grid_size, grid_size,2);
S_n1 = nan(grid_size, grid_size, 2);

% x in rows
% y in columns 
flow = repmat(grid,grid_size,1).*repmat(grid',1,grid_size)-b;
U = sum(u_n); 
V = phi*sum(v_n(:,:,1),'all') + (1-phi)*sum(v_n(:,:,2),'all');

S_n = S_init;

check_S=1;
it = 0;
while (check_S > tol)
    if it == MaxIt
        disp('surplus failed to converge!')
        break
    end 
    % define S+ = max(S(x,y,sig),0):
    S_plus = max(S_n,0);
    % define total number of unemployed and vacancies:
    S_n1(:,:,1) = (flow + (1-sigL)*beta*S_n(:,:,1) - (1-alpha)*beta*lambda*(1/U)*u_n*S_plus(:,:,1)...
        - alpha*beta*lambda*( phi*S_plus(:,:,1)*v_n(:,:,1)'*(1/V) ...
                    + (1-phi)*S_plus(:,:,2)*v_n(:,:,2)'*(1/V) ));
    S_n1(:,:,2) = (flow + (1-sigL)*beta*S_n(:,:,2) - (1-alpha)*beta*lambda*(1/U)*u_n*S_plus(:,:,2)...
        - alpha*beta*lambda*( phi*S_plus(:,:,1)*v_n(:,:,1)'*(1/V) ...
                    + (1-phi)*S_plus(:,:,2)*v_n(:,:,2)'*(1/V) ));

    check_S = max(abs(S_n1 - S_n),[],'all');
    
    S_n = S_n1;
    it = it + 1;
    
end

if check_S < tol
    disp('surplus converged')
end

end

