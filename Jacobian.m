function J = Jacobian(x)
    n = max(size(x));
    k = eye(n);
    %for j = 1:n
    for j = 1:n
    for i=1:n
        xad(i,j) = valder(x(i), k(i,:));
    end
    end
    
    y = func3d(xad(:,1));
    J = y.der;

    