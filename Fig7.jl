using ApproxFun
using ApproxFunBase
using LinearAlgebra
using GenericLinearAlgebra
using SpecialFunctions
using Plots
using LaTeXStrings

function example6(n , e)
    ## test for Airy equation
    ## eps * u''(x) - x * u(x) = 0
    ## u(-1) = Ai(-sqrt3(1 / eps))
    ## u(1) = Ai(sqrt3(1 / eps))
    ## u_exact = Ai(sqrt3(1 / eps) * x);
    ## solution

    x = Fun();
    eps = e;

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

    ## right hand side
    f = zeros(n , 1);
    u0 = S * u0;
    f[3 : end] = u0[1 : end - 2];
    f[1] = airyai(-1 / (eps)^(1/3));
    f[2] = airyai(1 / eps^(1/3));
    
    ## solve 
    # unum = L \ f;

    # Q , R = householderQR(L, 5)
    # unum = back_substitution(R, Q' * f)

    ~, R = householderQR([L f], 5);
    unum = back_substitution(R[: , 1 : end - 1], R[: , end])

    # Q , R = givensQR(L, 5)
    # unum = back_substitution(R, Q' * f) 

    return unum
end

function example8(n , e)
    ## test for Airy equation
    ## eps * u''(x) - x * u(x) = 0
    ## u(-1) = Ai(-sqrt3(1 / eps))
    ## u(1) = Ai(sqrt3(1 / eps))
    ## u_exact = Ai(sqrt3(1 / eps) * x);
    ## solution

    x = Fun();
    eps = e;

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

    ## right hand side
    f = zeros(n , 1);
    u0 = S * u0;
    f[3 : end] = u0[1 : end - 2];
    f[1] = airyai(-1 / (eps)^(1/3));
    f[2] = airyai(1 / eps^(1/3));
    
    ## solve 
    # unum = L \ f;

    # Q , R = householderQR(L, 5)
    # unum = back_substitution(R, Q' * f)

    Q, R = householderQR([L f], 5);
    unum = back_substitution(R[: , 1 : end - 1], R[: , end])
    # Q , R = givensQR(L, 5)
    # unum = back_substitution(R, Q' * f) 

    ## err
    uexact = coefficients(Fun(airyai(1 / eps^(1/3) * x) , Chebyshev() , n));
    err = norm(unum - uexact , 2);  ## solution(Float64) - solution(exact)
    # F = qr(L);
    # unum = F.R\(F.Q \ f);
    # unum = F \ f;

    # err
    return norm(Q[1 , end - 3 : end] , 2) + norm(Q[2 , end - 3 : end] , 2)
end

function back_substitution(U, b)
    n = size(U, 1)
    x = zeros(n)
    for i = n:-1:1
        if U[i, i] == 0
            error("矩阵是奇异的，不能求解")
        end
        x[i] = b[i]
        for j = i+1:n
            x[i] -= U[i, j] * x[j]
        end
        x[i] /= U[i, i]
    end
    return x
end

function householderQR(A, k)
    m, n = size(A)
    Q = Matrix{Float64}(I, m, m)  
    R = copy(A)  

    for i = 1:m
        x = R[i:min(m,i+k), i]

        e1 = zeros(length(x))
        e1[1] = 1
        v = sign(x[1]) * norm(x) * e1 + x
        v /= norm(v)  

        R[i:min(m,i+k), i:n] -= 2 * v * (v' * R[i:min(m,i+k), i:n])

        Q[:, i:min(m,i+k)] -= 2 * (Q[:, i:min(m,i+k)] * v) * v'
    end

    R = triu(R)
    return Q, R
end

function givensQR(A, k)
    m, n = size(A)
    Q = Matrix{Float64}(I, m, m)

    for i = 1:n-1
        for j = i+1:min(m, i+k)
            x = A[i:min(m, i+k), i]
            rt = givens(x, j-i+1)
            Q[:, i:min(m, i+k)] = Q[:, i:min(m, i+k)] * rt'
            A[i:min(m, i+k), :] = rt * A[i:min(m, i+k), :]
            # println("j: ", j)
        end
    end

    return Q, A
end

function givens(x, j)
    xi = x[1]
    xj = x[j]
    r = sqrt(xi^2 + xj^2)
    cost = xi / r
    sint = xj / r
    R = Matrix{Float64}(I, length(x), length(x))
    R[1, 1] = cost
    R[1, j] = sint
    R[j, 1] = -sint
    R[j, j] = cost
    return R
end

# Plots
x = 10 : 10 : 1000;
lx = size(x);
e1 = 1e-2;
global errcauchy1 = zeros(lx[1] , 1);
global err1 = zeros(lx[1] , 1);
for i = 1 : 1 : lx[1]
    u1 = example6(x[i], e1);
    u2 = example6(Int(ceil(1.01* x[i])), e1);
    global errcauchy1[i] = norm(u1 - u2[1 : x[i]] , 2);
    global err1[i] = example8(x[i], e1);
end
p1 = plot(x, log10.(errcauchy1), color=:black, shape=:circle, linewidth=2, xlabel=L"n", label=:"Cauchy error 1e-2", alpha=1, grid=false)
p2 = plot!(x, log10.(err1), color=:black, shape=:star4, linewidth=2, xlabel=L"n", label=:"Q 1e-2", alpha=1, grid=false)

e2 = 1e-5;
global errcauchy2 = zeros(lx[1] , 1);
global err2 = zeros(lx[1] , 1);
for i = 1 : 1 : lx[1]
    u1 = example6(x[i], e2);
    u2 = example6(Int(ceil(1.01* x[i])), e2);
    global errcauchy2[i] = norm(u1 - u2[1 : x[i]] , 2);
    global err2[i] = example8(x[i], e2);
end

p3 = plot!(x, log10.(errcauchy2), color=:black, shape=:star5, linewidth=2, xlabel=L"n", label=:"Cauchy error 1e-5", alpha=1, grid=false)
p4 = plot!(x, log10.(err2), color=:black, shape=:star6, linewidth=2, xlabel=L"n", label=:"Q1e-5", alpha=1, grid=false)


open("Fig7.txt", "w") do file
    for i in 1:length(x)
        write(file, "$(errcauchy1[i])\t$(err1[i])\t$(errcauchy2[i])\t$(err2[i])\t$(x[i])\n")
    end
end
