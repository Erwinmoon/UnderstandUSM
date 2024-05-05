using ApproxFun
using ApproxFunBase
using LinearAlgebra
using GenericLinearAlgebra
using SpecialFunctions
using Plots
using LaTeXStrings

function example10(n)
    ## test for Airy equation
    ## eps * u''(x) - x * u(x) = 0
    ## u(-1) = Ai(-sqrt3(1 / eps))
    ## u(1) = Ai(sqrt3(1 / eps))
    ## u_exact = Ai(sqrt3(1 / eps) * x);
    ## solution

    x = Fun();
    eps = 1e-5;

    ## innital date u0bound
    u0 = zeros(n , 1);

    ## diff mat D
    D = zeros(n , n);
    D[1:n, 1:n] = Derivative(Chebyshev() , 2)[1:n, 1:n];

    ## convert matrix S
    S = zeros(n , n);
    S = Conversion(Chebyshev(),Ultraspherical(2))[1:n , 1:n];

    ## multimat M
    M = zeros(n , n);
    M = Multiplication(x, Chebyshev())[1:n , 1:n];

    ## boundary conditions B
    B = zeros(2 , n);
    B = Dirichlet(-1..1)[1 : 2 , 1 : n];

    ## diff operator
    L = zeros(n , n);
    L = eps .* D  - S * M;
    L[3 : end , :] = L[1 : end - 2 , :];
    L[1 : 2 , :] = B;
    
    ## solve 
    # Q, R = givensQR(L, 5);
    F = qr(L);
    return norm(F.Q[1:2,n-3:n] , 2)
end

# Plots
x = 10 : 10 : 350;
lx = size(x);
global errcauchy = zeros(lx[1] , 1);
global err = zeros(lx[1] , 1);
for i = 1 : 1 : lx[1]
    global err[i] = example10(x[i]);
end

p1 = plot(x, log10.(err), color=:black, shape=:circle, linewidth=2, xlabel=L"n", label=:"Cauchy error", alpha=1, grid=false)