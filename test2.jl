using ApproxFun
using ApproxFunBase
using LinearAlgebra
using GenericLinearAlgebra
using SpecialFunctions
## using Plots

## some command
## new command : show all elements of A(array) --> show(stdout, "text/plain", A)
## plot(x , log10.(err), color=:blue, linewidth=2, xlabel="n", ylabel="err",label="", title="Float64-BigFloat",yticks=[-17,-16,-15,-14], alpha=0)

function example1(n)
    ## test for Airy equation
    ## eps * u''(x) - x * u(x) = 0
    ## u(-1) = Ai(-sqrt3(1 / eps))
    ## u(1) = Ai(sqrt3(1 / eps))
    ## u_exact = Ai(sqrt3(1 / eps) * x);
    ## cauchy error

    x = Fun();
    ep = 1e-2;
    p = Inf; ## p-norm

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
    unum = qr(L) \ f;
    bunum = bL \ bf;
    uexact = coefficients(Fun(airyai(1 / ep^(1/3) * x) , Chebyshev() , n));
    err1 = norm(unum - uexact , p);  ## solution(Float64) - solution(exact)
    err2 = norm(bunum - uexact , p);  ## solution(BigFloat) - solution(exact)
    err3 = norm(bunum - unum , p);  ## solution(Float64) - solution(BigFloat)

    ## defind m, such that f[m + 1 ; end] < e , (unit rounding error)
    e = 2^(-52);
    global m = 0;
    for ii = 1 : 1 : n
        if f[ii] > e
            global m = ii;
        end
    end
    invL = inv(L);
    AA = abs.(invL) * abs.(L);
    # bound = e * norm(invL[: , 1 : m] , p) * norm(f[1 : m] , p) + e^2 * norm(invL[: , m + 1 : end] , p) * norm(ones(n , 1)); ## norm wise boundary
    ker = abs.(invL[: , 1 : m]) * abs.(f[1 : m]) + AA[: , 1 : m] * abs.(bunum[1 : m])
    bound = e * norm(ker,p);
    ## conditon numbers
    condeff = CondEff(nL , unum , f) * norm(unum , p);
    condnorm = cond(L , p);

    # return err1, err2, err3, L , m , bound , f , unum
    return err3, err1, condeff, condnorm, bound
    # return unum
end

function example2(n)        
    ## test for 
    ## u'(x) = cos(x)
    ## u(-1) = sin(-1)
    ## u_exact = cos(x);

    x = Fun();

    ## innital date u0bound
    u0 = coefficients(Fun(cos(x) , Chebyshev() , n));;

    ## diff mat D
    D = zeros(n , n);
    D[1:n, 1:n] = Derivative(Chebyshev() , 1)[1:n, 1:n];


    ## convert matrix S
    S = zeros(n , n);
    S = Conversion(Chebyshev(),Ultraspherical(1))[1:n , 1:n];

    ## multimat M
    M = zeros(n , n);
    M = Multiplication(1, Chebyshev())[1:n , 1:n];

    ## boundary conditions B
    B = zeros(1 , n);
    B = Dirichlet(-1..1)[1 , 1 : n];

    ## diff operator
    L = zeros(n , n);
    L = D;
    L[2 : end , :] = L[1 : end - 1 , :];
    L[1 , :] = B;

    ## right hand side
    f = zeros(n , 1);
    u0 = S * u0;
    f[2 : end] = u0[1 : end - 1];
    f[1] = sin(-1);

    ## defind m, such that f[m + 1 ; end] < e , (unit rounding error)
    p = Inf; ## p-norm
    e = 2^(-52);
    global m = 0;
    for ii = 1 : 1 : n
        if f[ii] > e
            global m = ii;
        end
    end
    invL = inv(L);
    bound = 2 * e * norm(invL[: , 1 : m] , p) * norm(f[1 : 3] , p); ## norm wise boundary

    ## BigFloat 
    bL = convert(Matrix{BigFloat}, L);
    bf = convert(Matrix{BigFloat}, f);

    ## solve
    unum = L \ f;
    bunum = bL \ bf;
    uexact = coefficients(Fun(sin(x) , Chebyshev() , n));
    err1 = norm(unum - uexact , p);  ## solution(Float64) - solution(exact)
    err2 = norm(bunum - uexact , p);  ## solution(BigFloat) - solution(exact)
    err3 = norm(bunum - unum , p);  ## solution(Float64) - solution(BigFloat)

    return err1, err2, err3, L, m, bound, f    
end

function example3(n)
    ## test for 
    ## u'(x) + x^3 * u(x) = 100sin(20000x^2)
    ## u(-1) = 0
    ## u_exact = ?;
    ## n >= 315; 

    x = Fun();

    ## innital date u0bound
    u0 = coefficients(Fun(100 * sin(200 * x^2) , Chebyshev()));
    mm = size(u0);
    if mm[1] < n
        u0 = [u0 ; zeros(n - mm[1] , 1)];
    end
    ## diff mat D
    D = zeros(n , n);
    D[1:n, 1:n] = Derivative(Chebyshev() , 1)[1:n, 1:n];


    ## convert matrix S
    S = zeros(n , n);
    S = Conversion(Chebyshev(),Ultraspherical(1))[1:n , 1:n];

    ## multimat M
    M = zeros(n , n);
    M = Multiplication(x^3, Chebyshev())[1:n , 1:n];

    ## boundary conditions B
    B = zeros(1 , n);
    B = Dirichlet(-1..1)[1 , 1 : n];

    ## diff operator
    L = zeros(n , n);
    L = D;
    L[2 : end , :] = L[1 : end - 1 , :];
    L[1 , :] = B;

    ## right hand side
    f = zeros(n , 1);
    u = S * u0[1 : n];
    f[2 : end] = u[1 : n - 1];
    f[1] = 0;

    ## defind m, such that f[m + 1 ; end] < e , (unit rounding error)
    p = Inf; ## p-norm
    e = 2^(-51);
    global m = 0;
    for ii = 1 : 1 : n
        if f[ii] > e
            global m = ii;
        end
    end
    invL = inv(L);
    bound = 2 * e * norm(invL[: , 1 : m] , p) * norm(f[1 : 3] , p); ## norm wise boundary

    ## BigFloat 
    bL = convert(Matrix{BigFloat}, L);
    bf = convert(Matrix{BigFloat}, f);

    ## solve
    unum = L \ f;
    bunum = bL \ bf;
    uexact = bunum;
    err1 = norm(unum - uexact , p);  ## solution(Float64) - solution(exact)
    err2 = norm(bunum - uexact , p);  ## solution(BigFloat) - solution(exact)
    err3 = norm(bunum - unum , p);  ## solution(Float64) - solution(BigFloat)

    return err1, err2, err3, L, m, bound, f
end

function example4(n)
    ## test for Airy equation
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

    ## right hand side
    f = zeros(n , 1);
    u0 = S * u0;
    f[3 : end] = u0[1 : end - 2];
    f[1] = airyai(-1 / (eps)^(1/3));
    f[2] = airyai(1 / eps^(1/3));

    ## defind m, such that f[m + 1 ; end] < e , (unit rounding error)
    p = Inf; ## p-norm
    e = 2^(-52);
    global m = 0;
    for ii = 1 : 1 : n
        if f[ii] > e
            global m = ii;
        end
    end


    
    ## solve 
    unum = L \ f;
    uexact = coefficients(Fun(airyai(1 / eps^(1/3) * x) , Chebyshev() , n));
    err1 = norm(unum - uexact , p);  ## solution(Float64) - solution(exact)
    F = qr(L);
    # unum = F.R\(F.Q \ f);
    unum = F \ f;

    # return err1, L, unum, m, uexact
    return unum
end

function example5(n)        
    ## test for 
    ## u'(x)  + (1 / (ax^2 + 1))u = 0
    ## u(-1) = 1
    ## u_exact = exp(-[acrtan(sqrt(a)x) + arctan(sqrt(a))] / (sqrt(a)));

    a = 5e4;
    x = Fun();

    ## innital date u0bound
    u0 = coefficients(Fun(0 * x , Chebyshev() , n));

    ## diff mat D
    D = zeros(n , n);
    D[1:n, 1:n] = Derivative(Chebyshev() , 1)[1:n, 1:n];


    ## convert matrix S
    S = zeros(n , n);
    S = Conversion(Chebyshev(),Ultraspherical(1))[1:n , 1:n];

    ## multimat M
    M = zeros(n , n);
    M = Multiplication(1 / (a * x^2 + 1), Chebyshev())[1:n , 1:n];

    ## boundary conditions B
    B = zeros(1 , n);
    B = Dirichlet(-1..1)[1 , 1 : n];

    ## diff operator
    L = zeros(n , n);
    L = D + S * M;
    L[2 : end , :] = L[1 : end - 1 , :];
    L[1 , 1 : n] = B;

    ## right hand side
    f = zeros(n , 1);
    u0 = S * u0;
    f[2 : end] = u0[1 : end - 1];
    f[1] = 1;

    ## solve
    unum = L \ f;

    return unum, L   
end

function example6(n)
    ## test for Airy equation with a(x) = 1
    ## eps * u''(x) - u(x) = 0
    ## cauchy error

    x = Fun();
    eps = 1e-7;

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
    unum = L \ f;
    ## p = Inf;
    ## bL = convert(Matrix{BigFloat}, L);
    ## bf = convert(Matrix{BigFloat}, f);
    ## bunum = convert(Matrix{BigFloat}, unum);
    ## niu = norm(bL * bunum - bf , p) / (norm(L , p) * norm(unum , p) + norm(f , p)); ## normwise backward 
    ## omi = norm((bL * bunum - bf) ./ (abs.(L) * abs.(unum) + f) , Inf);
    ## res = norm(bL * bunum - bf , p);
    return unum
end

function example7(n)
    ## test for Airy equation with a(x) = 1
    ## eps * u''(x) - a(x) * u'(x) - b(x) * u(x) = 0
    ## cauchy error

    x = Fun();
    eps = 1e-0;

    ## inital date u0
    u0 = zeros(n , 1);

    ## diff mat D
    D2 = zeros(n , n);
    D2[1:n, 1:n] = Derivative(Chebyshev() , 2)[1:n, 1:n];
    D1 = zeros(n , n);
    D1[1:n, 1:n] = Derivative(Chebyshev() , 2)[1:n, 1:n];

    ## convert matrix S
    Sb = zeros(n , n);
    Sb = Conversion(Chebyshev(),Ultraspherical(2))[1:n , 1:n];
    Sa = zeros(n , n);
    Sa = Conversion(Ultraspherical(1),Ultraspherical(2))[1:n , 1:n];

    ## multimat M
    Mb = zeros(n , n);
    ## Mb = Multiplication(0 * x, Chebyshev())[1:n , 1:n];  ## M[b(x)]
    Ma = zeros(n , n);
    Ma = Multiplication(1 * x, Ultraspherical(1))[1:n , 1:n];  ## M[a(x)]

    ## boundary conditions B
    B = zeros(2 , n);
    B = Dirichlet(-1..1)[1 : 2 , 1 : n];

    ## diff operator
    L = zeros(n , n);
    L = eps .* D2  + Sa * D1 * Ma - Sb * Mb;
    L[3 : end , :] = L[1 : end - 2 , :];
    L[1 : 2 , :] = B;

    ## right hand side
    f = zeros(n , 1);
    u0 = Sb * u0;
    f[3 : end] = u0[1 : end - 2];
    f[1] = airyai(-1 / (eps)^(1/3));
    f[2] = airyai(1 / eps^(1/3));


    ## solve
    unum = L \ f;
    ## p = Inf;
    ## bL = convert(Matrix{BigFloat}, L);
    ## bf = convert(Matrix{BigFloat}, f);
    ## bunum = convert(Matrix{BigFloat}, unum);
    ## niu = norm(bL * bunum - bf , p) / (norm(L , p) * norm(unum , p) + norm(f , p)); ## normwise backward 
    ## omi = norm((bL * bunum - bf) ./ (abs.(L) * abs.(unum) + f) , Inf);
    ## res = norm(bL * bunum - bf , p);
    return unum
end

function example8(n)
    ## example of 10-order ode
    x = Fun();
    p = 4;
    order = 2 * (p + 1);

    ## some special matrices
    D10 = zeros(n , n);
    D10[1:n, 1:n] = Derivative(Chebyshev() , order)[1:n, 1:n];
    S10 = zeros(n , n);
    S10 = Conversion(Chebyshev(),Ultraspherical(order))[1:n , 1:n];

    l = -1;
    r = 1;
    ## boundary conditons
    global B = zeros(2 * (p + 1) , n);
    for i = 1 : 1 : p + 1
       global B[2 * i - 1 , :] = HighOrderBC(n , i - 1 , l);
    end
    for i = 1 : 1 : p + 1
       global B[2 * i , :] = HighOrderBC(n , i - 1 , r);
    end

    # LHS
    L = zeros(n , n);
    L[1 : 2 * (p + 1) , :] = B;
    L[2 * (p + 1) + 1 : end , :] = D10[1 : n - 2 * (p + 1) , :];

    # RHS
    f = zeros(n , 1);
    f = S10 * coefficients(Fun(x -> -sin(x) , Chebyshev() , n));
    f[2 * (p + 1) + 1 : end] = f[1 : n - 2 * (p + 1)];
    f[1] = sin(l);
    f[2] = sin(r);
    f[3] = cos(l);
    f[4] = cos(r);
    f[5] = -sin(l);
    f[6] = -sin(r);
    f[7] = -cos(l);
    f[8] = -cos(r);
    f[9] = sin(l);
    f[10] = sin(r);

    # solve
    bL = convert(Matrix{BigFloat}, L);
    bf = convert(Vector{BigFloat}, f);
    u = qr(L) \ f;
    bu = L\bf;
    uexact = coefficients(Fun(x -> sin(x) , Chebyshev() , n));
    err = norm(u - bu , Inf);
    err1 = norm(u - uexact , Inf);
    condeff = CondEff(L , u , f) * norm(bu , Inf);
    condnorm = cond(L , Inf);

    e = 2^(-51);
    global m = 0;
    for ii = 1 : 1 : n
        if f[ii] > e
            global m = ii;
        end
    end
    invL = inv(L);
    bound = 2 * e * norm(invL[: , 1 : m] , p) * norm(f[1 : m] , p) + 2 * e^2 * norm(invL[: , m + 1 : end] , p) * norm(ones(n , 1));

    # return err, err1, condeff, condnorm, bound
    return u;
end

function example9(n)
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
    c = cond(L , 2);

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

function CondEff(A , x , b)
    p = 2;
    F = svd(A);
    c = norm(A*x - b , p) / (F.S[end] * norm(x , p));
    return c
end

# main function
 # n = 1000;
 # (err1, err2, err3, L , m , bound , f , unum) = example1(n);
 # unum = example4(n);



    ## main fuction : test for cauchy error for airy equation
    #(~, L1, unum1, ~) = example4(n);
    #m = Int(ceil(1.01 * n));
    #(err1, L2, unum2, ~ , uexact) = example4(m);
    #err = norm(unum1 - unum2[1 : n] , 2)
    #bound = norm(L[1 : n , 1 : n] \ L[1 : n , n + 1 : m] * unum2[n + 1 : m] , 2);
#=
    ## main function : polt err1
    using Plots
    x = 10 : 10 : 500;
    lx = size(x);
    global err = zeros(lx[1] , 1);
    global bound = zeros(lx[1] , 1);
    for i = 1 : 1 : lx[1]
        n = x[i];
        (~, ~, err1, ~ , ~ , bound1 , ~ , ~) = example1(n);
        global err[i] = err1;
        global bound[i] = bound1;
    end
    plot(x , [err bound], shape = [:circle :cross],color=:black, linewidth=2, xlabel="n", ylabel="condition number",label=["error" "bound"], alpha=1)
    # savefig(p, "name.pdf")  ## save fig
    # plot!(legend=:outerbottom, legendcolumns=3)  \\  plot!(legend=:right) \\ grid=false 

    ## err and bound 
    using Plots
    p = Inf; ## p-norm
    x = 100 : 10 : 500;
    lx = size(x);
    global err = zeros(lx[1] , 1);
    global err1 = zeros(lx[1] , 1);
    global err2 = zeros(lx[1] , 1);
    global err3 = zeros(lx[1] , 1);
    global err4 = zeros(lx[1] , 1);
    e = 2^(-51); 
    for i = 1 : 1 : lx[1]
        local ans = example1(x[i]);
        global err[i] = ans[1];
        global err1[i] = ans[2]; 
        global err2[i] = ans[3]; 
        global err3[i] = e * ans[4]; 
        global err4[i] = ans[5];
    end
    # plot(x, [log10.(err),log10.(err1)], shape = [:circle :cross],color=:black, linewidth=2, xlabel="n", ylabel= "err" ,label=["Float64-BigFloat" "Float64-uexact"], alpha=1)
    # plot(x, [log10.(err),log10.(err1), log10.(err2), log10.(err3), log10.(err4)], shape = [:circle :cross :star5 :star4 :star2],color=:black, linewidth=2, xlabel="n", ylabel= "err" ,label=["Float64-BigFloat" "Float64-uexact" "condeff" "cond" "mycond"], alpha=1)
    # plot(x, [err,err1, err2, err4], shape = [:circle :cross :star5 :star4 ],color=:black, linewidth=2, xlabel="n", ylabel= "err" ,label=["Float64-BigFloat" "Float64-uexact" "condeff" "mycond"], alpha=1)
 
         ## err and bound 
    using Plots
    p = Inf; ## p-norm
    x = 100 : 10 : 200;
    lx = size(x);
    global err = zeros(lx[1] , 1);
    global err1 = zeros(lx[1] , 1);
    global err2 = zeros(lx[1] , 1);
    global err3 = zeros(lx[1] , 1);
    global err4 = zeros(lx[1] , 1);
    e = 2^(-51); 
    for i = 1 : 1 : lx[1]
        local ans = example1(x[i]);
        global err[i] = ans[1];
        global err1[i] = ans[2]; 
        global err2[i] = ans[3]; 
        global err3[i] = e * ans[4]; 
        global err4[i] = ans[5];
    end
    # plot(x, [log10.(err),log10.(err1)], shape = [:circle :cross],color=:black, linewidth=2, xlabel="n", ylabel= "err" ,label=["Float64-BigFloat" "Float64-uexact"], alpha=1)
    plot(x, [log10.(err),log10.(err1), log10.(err2), log10.(err3), log10.(err4)], shape = [:circle :cross :star5 :star4 :star2],color=:black, linewidth=2, xlabel="n", ylabel= "err" ,label=["Float64-BigFloat" "Float64-uexact" "condeff" "cond" "mycond"], alpha=1)
    # plot(x, [err,err1, err2, err4], shape = [:circle :cross :star5 :star4 ],color=:black, linewidth=2, xlabel="n", ylabel= "err" ,label=["Float64-BigFloat" "Float64-uexact" "condeff" "mycond"], alpha=1)

    

    ## main function : example in "A Fast ..." Fig. 2.2
    using Plots
    p = Inf; ## p-norm
    x = 100 : 10 : 500;
    lx = size(x);
    global err = zeros(lx[1] , 1);
    for i = 1 : 1 : lx[1]
        u1 = example4(x[i]);
        u2 = example4(Int(ceil(1.01* x[i])));
        global err[i] = norm(u1 - u2[1 : x[i]] , p);
    end
    plot(x , log10.(err), color=:black, linewidth=2, xlabel="n", ylabel="cauchy error",label=:none,alpha=1)
    =#


    using Plots
    p = Inf; ## p-norm
    x = 100 : 10 : 500;
    lx = size(x);
    global err = zeros(lx[1] , 1);
    for i = 1 : 1 : lx[1]
        global err[i] = example9(x[i]);
    end
    plot(x , err, color=:black, linewidth=2, xlabel="n", ylabel="conditon number",label=:none,alpha=1)



    