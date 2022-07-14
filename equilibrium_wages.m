function w = equilibrium_wages(S,u_n,v_n,b,alpha,beta,sigL,sigH,phi,lambda,grid_size,tol,tol_out,MaxIt)

w = nan(grid_size,grid_size,2);
S_plus = max(S,0);
V = phi*sum(v_n(:,:,1),'all') + (1-phi)*sum(v_n(:,:,2),'all');

for x = 1:grid_size
   for y = 1:grid_size
       w(x,y,1) = (1- beta*(1-sigL))*alpha*S(x,y,1) + b ...
                   + alpha*beta*lambda*( phi*S_plus(x,:,1)*v_n(:,:,1)'/V ...
                        + (1-phi)*S_plus(x,:,2)*v_n(:,:,2)'/V) ;
       w(x,y,2) = (1- beta*(1-sigH))*alpha*S(x,y,2) + b ...
                   + alpha*beta*lambda*( phi*S_plus(x,:,1)*v_n(:,:,1)'/V ...
                        + (1-phi)*S_plus(x,:,2)*v_n(:,:,2)'/V) ;
   end
end 

end