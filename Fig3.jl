using ApproxFun
using ApproxFunBase
using LinearAlgebra
using GenericLinearAlgebra
using SpecialFunctions
using Plots
using LaTeXStrings

function example7(n)
    ## test for Airy equation
    ## eps * u''(x) - x * u(x) = 0
    ## u(-1) = Ai(-sqrt3(1 / eps))
    ## u(1) = Ai(sqrt3(1 / eps))
    ## u_exact = Ai(sqrt3(1 / eps) * x);
    ## cauchy error

    x = Fun();
    ep = 1e-2;
    p = Inf; ## p-norm
    m = 3;

    ## inital date u0
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
    l = -1;
    r = 1;
    global B = zeros(2  , n);
    for i = 1 : 1 : 1
       global B[2 * i - 1 , :] = HighOrderBC(n , i - 1 , l);
    end
    for i = 1 : 1 : 1
       global B[2 * i , :] = HighOrderBC(n , i - 1 , r);
    end

    ## diff operator
    L = zeros(n , n);
    L = ep .* D  - S * M;
    L[3 : end , :] = L[1 : end - 2 , :];
    L[1 : 2 , :] = B;

    ## right hand side
    f = zeros(n , 1);
    u0 = S * u0;
    f[3 : end] = u0[1 : end - 2];
    f[1] = airyai(-1 / (ep)^(1/3));
    f[2] = airyai(1 / ep^(1/3));

    ## BigFloat 
    bL = convert(Matrix{BigFloat}, L);
    bf = convert(Matrix{BigFloat}, f);
    nL = convert(Matrix{Float64}, L);

    ## solve 
    # unum = L \ f;

    Q , R = householderQR(L, 5)
    unum = back_substitution(R, Q' * f)

    # Q , R = givensQR(L, 5)
    # unum = back_substitution(R, Q' * f) 
    
    bunum = bL \ bf;

    err = norm(bunum - unum , p);  ## solution(Float64) - solution(BigFloat)

    ## defind m, such that f[m + 1 ; end] < e , (unit rounding error)
    e = 2^(-52);
    global r = 0;
    for ii = 1 : 1 : n
        if unum[ii] > e
            global r = ii;
        end
    end 
    invL = inv(L);
    AA = abs.(invL[ : , 1 : min(r + m, n)]) * abs.(L[1 : min(r + m, n) , 1 : r]);
    # mybound = e * norm(invL[: , 1 : m] , p) * norm(f[1 : m] , p) + e^2 * norm(invL[: , m + 1 : end] , p) * norm(ones(n , 1)); ## norm wise boundary
    ker = abs.(invL[: , 1 : r]) * abs.(f[1 : r]) + AA[: , 1 : r] * abs.(bunum[1 : r]);
    mybound = e * norm(ker , p) / norm(unum , p);
    ## conditon number bound 
    condbound = e * cond(L , p);

    return err, mybound, condbound
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

p = 2; ## p-norm
x = 2 : 2 : 350;
lx = size(x);
global err = zeros(lx[1] , 1);
global mybound = zeros(lx[1] , 1);
global cbound = zeros(lx[1] , 1);
for i = 1 : 1 : lx[1]
    local ans = example7(x[i]);
    global err[i] = ans[1];
    global mybound[i] = ans[2]; 
    global cbound[i] = ans[3]; 
end
plot(x, [log10.(err),log10.(mybound), log10.(cbound)], shape = [:circle :cross :star5],color=:black, linewidth=2, xlabel=L"n",label=["Float64-BigFloat" "mycond" "cond"], alpha=1, grid=false)
plot!(legend=:right)

# open("Fig3.txt", "w") do file
   #  for i in 1:length(x)
   #      write(file, "$(err[i])\t$(mybound[i])\t$(cbound[i])\t$(x[i])\n")
   #  end
# end