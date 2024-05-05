clear
clc


nn = 10 : 10 : 100;
err = zeros(1 , length(nn));

for i = 1 : 1 : length(nn)
    err(i) = Example(nn(i));
end
semilogy(nn, err)

function err = Example(n)
    N = chebop2(@(u) diff(u,2,1) + diff(u,2,2));
    f = chebfun2(@(x,y) 2*x^2 + 2*y^2 - 4, 'vectorize'); 
    N.bc = 0;
    u = chebop2.denseSolve(N, f, n, n);
    
    uexact = chebfun2(@(x,y) (x^2 - 1) * (y^2 - 1), 'vectorize');
    err = norm(u - uexact,inf);
end