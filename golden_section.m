function x_min = golden_section(fun,x1,x2)

%x1 is the start of interval; x2 is the end of interval
Tolx= (x2-x1)*1e-3; %Tolerance
iter = 0; %Number of iterations
maxiter = 100; %Maximal iterations
t = double((sqrt(5)-1)/2); %Interval reduced factor

x1_new = x1 + (1-t)*(x2-x1);
x2_new = x1 + t*(x2-x1);
f_1 = fun(x1_new);
f_2 = fun(x2_new);

fprintf('  iteration number   x    f(x)\n')
data= [];
while (abs(x2-x1)>Tolx) && (iter<maxiter)
 
 if(f_1<f_2)
     x2 = x2_new;
     x2_new = x1_new;
     x1_new = x1 + (1-t)*(x2-x1);
     f_1 = fun(x1_new);
     f_2 = fun(x2_new);
 else
     x1 = x1_new;
     x1_new = x2_new;
     x2_new = x1 + t*(x2-x1);
     f_1 = fun(x1_new);
     f_2 = fun(x2_new);
 end
    iter = iter+1;
    x_min = min(x1_new,x2_new);
    f_min = min(f_1,f_2);
    data = [data;iter,x_min,f_min];
    
end
disp(data)

