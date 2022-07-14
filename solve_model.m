function [S_init,u_n,v_n] = solve_model(b,alpha,beta,sigL,sigH,phi,lambda,grid_size,tol,tol_out,MaxIt)

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
it_out = 0;
while (check > tol_out)
    if it_out == MaxIt
        disp('Failed to converge!')
        break
    end 
    it_out = it_out + 1

    S_n1 = VFI_surplus(S_init);

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

    
end; if check < tol_out; disp('converged.'); end 






function S_n = VFI_surplus(S)
% -------------------------------------------------------------------------
% This function performs a value function ...
% 
% Inputs:
% - S: initial guess for surplus function 
% - ... other parameters from exercise 
% - u: density of unemployed (1xgrid_size) vector
% - v: density of vacancies  (1xgrid_sizex2) vector
% - tol: tolerance value for convergence 
% - MaxIt: Maximum number of iterations 
% -------------------------------------------------------------------------

% Creating equidistant grid for productivities x,y
grid = linspace(0,1,grid_size);

S_plus = nan(grid_size, grid_size,2);
S_new = nan(grid_size, grid_size, 2);

% x in rows
% y in columns 
flow = repmat(grid,grid_size,1).*repmat(grid',1,grid_size)-b;
U = sum(u_n); 
V = phi*sum(v_n(:,:,1),'all') + (1-phi)*sum(v_n(:,:,2),'all');

S_n = S;

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
    S_new(:,:,1) = (flow + (1-sigL)*beta*S_n(:,:,1) - (1-alpha)*beta*lambda*(1/U)*u_n*S_plus(:,:,1)...
        - alpha*beta*lambda*( phi*S_plus(:,:,1)*v_n(:,:,1)'*(1/V) ...
                    + (1-phi)*S_plus(:,:,2)*v_n(:,:,2)'*(1/V) ));
    S_new(:,:,2) = (flow + (1-sigL)*beta*S_n(:,:,2) - (1-alpha)*beta*lambda*(1/U)*u_n*S_plus(:,:,2)...
        - alpha*beta*lambda*( phi*S_plus(:,:,1)*v_n(:,:,1)'*(1/V) ...
                    + (1-phi)*S_plus(:,:,2)*v_n(:,:,2)'*(1/V) ));

    check_S = max(abs(S_new - S_n),[],'all');
    
    S_n = S_new;
    it = it + 1;
    
end

if check_S < tol
    disp('surplus converged')
end

end

end