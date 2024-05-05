using ApproxFun
using ApproxFunBase
using LinearAlgebra
using GenericLinearAlgebra
using SpecialFunctions
using Plots
using LaTeXStrings

function example1(n)
    ## test for Airy equation condition number
    ## eps * u''(x) - x * u(x) = 0
    ## u(-1) = Ai(-sqrt3(1 / eps))
    ## u(1) = Ai(sqrt3(1 / eps))
    ## u_exact = Ai(sqrt3(1 / eps) * x);
    ## solution

    x = Fun();
    eps = 1e-2;

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

    ## cond
    c = cond(L , Inf);

    # return err1, L, unum, m, uexact
    return c
end

function example2(n)
    ## test for Airy equation condition number
    ## eps * u''(x) - x * u(x) = 0
    ## u(-1) = Ai(-sqrt3(1 / eps))
    ## u(1) = Ai(sqrt3(1 / eps))
    ## u_exact = Ai(sqrt3(1 / eps) * x);
    ## solution

    x = Fun();
    eps = 1e-2;

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
    B[1 , : ] = HighOrderBC(n , 0 , 1);
    B[2 , : ] = HighOrderBC(n , 1 , -1);

    ## diff operator
    L = zeros(n , n);
    L = eps .* D  - S * M;
    L[3 : end , :] = L[1 : end - 2 , :];
    L[1 : 2 , :] = B;

    ## cond
    c = cond(L , Inf);

    # return err1, L, unum, m, uexact
    return c
end

function HighOrderBC(n , p , point)
    B = zeros(1 , n);
    for i = 1 : 1 : n
        b = (point)^(i + p - 1);
        for k = 0 : 1 : p - 1
            b = b * ((i-1)^2 - k^2) / (2 * k + 1);
        end
        B[1 , i] = b;
    end
    return B
end

# Plots
x = 10 : 10 : 350;
lx = size(x);
global errdd = zeros(lx[1] , 1);
# global errdn = zeros(lx[1] , 1);
for i = 1 : 1 : lx[1]
    global errdd[i] = example1(x[i]);
    # global errdn[i] = example2(x[i]);
end
# 1st fig
# p1 = plot(x, errdd, color=:black, linewidth=2, xlabel=L"n", ylabel="condition number", label=:none, alpha=1, grid=false)

# 2ed fig
# cp2 = plot(x, [errdd, errdn], shape=[:circle :cross], color=:black, linewidth=2, xlabel=L"n", ylabel="condition number", label=["Dirichlet-Dirichlet" "Dirichlet-Neumann"], alpha=1, grid=false)

# layout = @layout [a; b]
# p = plot(p2, p1, layout=layout, size=(600, 400))

# display(p)

# open("Fig1a.txt", "w") do file
    # for i in 1:length(x)
        # write(file, "$(errdd[i])\t$(x[i])\n")
    # end
# end