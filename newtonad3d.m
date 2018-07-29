function [y, iter] = newtonad3d(x)
% [y, iter]= Newtonsystem(x) 
%Solve 3 by 3 system of nonlinear equations by Newton's method 
% Input x is intital guess to soln, a vector of dimension 3; 
% func is the vector function whose root is to be found, 
% Jacobian and func compute Jacobian and func values from input x 
% these should be defined in separate function files
% relative tolerance is hard-wired to 1e-8; max iterations to 20

x= x(:); % make a column vector
xold = x; xnew= x; iter = 0;  p = ones(3,1); % initial p is arbitrary

fprintf('  iteration number   x(1) x(2) x(3) f(1) f(2) f(3) p\n')
data= [];
while ((norm(p, 1) >= 1e-8 * (1 + norm(xold,1))) & iter <= 20) 
    
    
    F= func3d([valder(xnew(1),1);valder(xnew(2),1);valder(xnew(3),1)]); 
    F = F.val;
    J = Jacobian(xnew); 
    p= J\F; % solve for the step in Newton's method 
    xold = xnew; 
    xnew = xold - p; % Newton's method 
    iter = iter + 1;   
    data = [data;iter,xold(1),xold(2),xold(3),F(1),F(2),F(3),norm(p,1)];
    
end
disp(data)
y= xold
