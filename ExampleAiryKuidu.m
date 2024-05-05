clear
clc

%%
% base of CTM bu kuidu 
% example for some easy function 
% airy equation \eps * u'' - xu = 0

nn = 10 : 20 : 350;
errcauchy = zeros(1 , length(nn));
err = zeros(1 , length(nn));

for i = 1 : 1 : length(nn)
    n1 = nn(i);
    n2 = ceil(1.01 * n1);

    u1 = ExampleAiryCauchy(n1);
    u2 = ExampleAiryCauchy(n2);
    
    e = u1 - u2(1 : n1);

    errcauchy(i) = norm(e , 2);
    err(i) = ExampleAiryError(n1);
end

M = [nn' err' errcauchy'];
% writematrix(M, 'kuidu.txt', 'Delimiter', ',');
semilogy(nn, errcauchy)
hold on
semilogy(nn, err)

%%
% function for Cauchy error of airy equation
function u = ExampleAiryCauchy(n)
    e = 1e-2;
    m = 2; % order
    F = {chebfun(@(x)(-x)) , chebfun(0)};
    
    % LHS
    L = eye(n);
    Q = Qmatrix(n);
    A = Amatrix(n , m , F); 
    X = Xmatrix(n , m);
    
    
        % bc 
        B = zeros(m , n);
        B(1 , : ) = HighOrderBC(n , 0 , 1);
        B(2 , : ) = HighOrderBC(n , 0 , -1);
        b = zeros(m , 1);
        b(1) = airy(1 / (e^(1 / 3)));
        b(2) = airy(-1 / (e^(1 / 3)));
    
    M = ultraS.multmat(n , chebfun(@(x) -x) , 0);
    L = e * L + M * Q^m - A * inv(B * X) * B * Q^(m);
    
    % RHS
    f = zeros(n , 1);
    f = f - A * inv(B * X) * b;
    
    % solve
    [Q , R] = qr([L f]);
    % u = R(: , 1 : n) \ (R(: , end));
    u = upperTriangularSolve(R(: , 1 : n), R(: , end));
    % u = L \ f;
end

%%
% function for error airy equation ,u64 - uexact
function err = ExampleAiryError(n)
    e = 1e-2;
    m = 2; % order
    F = {chebfun(@(x)(-x)) , chebfun(0)};
    
    % LHS
    L = eye(n);
    Q = Qmatrix(n);
    A = Amatrix(n , m , F); 
    X = Xmatrix(n , m);
    
    
        % bc 
        B = zeros(m , n);
        B(1 , : ) = HighOrderBC(n , 0 , 1);
        B(2 , : ) = HighOrderBC(n , 0 , -1);
        b = zeros(m , 1);
        b(1) = airy(1 / (e^(1 / 3)));
        b(2) = airy(-1 / (e^(1 / 3)));
    
    M = ultraS.multmat(n , chebfun(@(x) -x) , 0);
    L = e * L + M * Q^m - A * inv(B * X) * B * Q^(m);
    
    % RHS
    f = zeros(n , 1);
    f = f - A * inv(B * X) * b;
    
    % solve
    % [Q , R] = qr(L);
    % v = R \ (Q' * f);
    % v = upperTriangularSolve(R, Q' * f);
    v = L \ f;
    
    u = Q^m * v + X * inv(B * X) * (b - B * Q^m * v);

    % error
    % D2 = ultraS.diffmat(n , 2);
    uexact = chebcoeffs(chebfun(@(x)airy(x / (e^(1 / 3))) , n)); % exact solution of e = 1 
    err = norm(u - uexact , inf);
end

%%
% test by D_k * Q^k = S_{k - 1} * ...* S_0
function value = test1(n)
    Q = Qmatrix(n);
    D2 = ultraS.diffmat(n , 2);
    S0 = ultraS.convertmat(n, 0 , 0);
    S1 = ultraS.convertmat(n, 1 , 1);
    A = D2 * Q * Q - S1 * S0;
    if norm(A(1 : ceil(n / 2) , 1 : ceil(n / 2)) , 2) <= eps
        value = 1;
    else 
        value = 0;
    end
end

%%
% test by eq between eq(1.3) and eq(1.4) in duiku
% run follow code
    % n = 10;
    % m = 2;
    % F = {chebfun(@(x)sin(x) , n), chebfun(@(x)cos(x) , n)};
    % v = test2(n , m, F);
function value = test2(n , m , F)
    x = chebfun(@(x)x , n);
    A = Amatrix(n , m , F);
    err1 = norm(chebcoeffs(F{1}) - A(: , 1) , "inf");
    err2 = norm(chebcoeffs(F{1} * x + F{2} , n) - A(: , 2) , "inf");
    if err1 < eps && err2 < eps
        value = 1;
    else 
        value = 0;
    end  
end

%%
% Qmatrix, cref[19] in dukui
% Q \in R^{n * n}
function Q = Qmatrix(n)
    Q = zeros(n , n);
    a = zeros(n , 1);
    for i = 2 : 1 : n - 1
        Q(1 , i + 1) = (-1)^(i + 1) / ((i - 1) * (i + 1));
    end
    for i = 1 : 1 : n
        a(i + 1) = 1 / (2 * i);
    end
    Q = Q + diag(a(2 : n) , -1);
    Q = Q + diag(-a(1 : n - 1) , 1);
    Q(1 , 1) = 1;
    Q(2 , 1) = 1;
    Q(1 , 2) = -1 / 4;
end

%%
% X(x), the matrix coefficients of [1 , x^1 , x^2, ...]
% X \in R^{n * m}
function X = Xmatrix(n , m)
    X = zeros(n , m);
    for i = 1 : 1 : m
        X(: , i) = chebcoeffs(chebfun(@(x) x^(i - 1)) , n);
    end
end

%%
% A(x), the matrix coefficients of [a^0 , a^0*x + a^1,...]
% A(x) \in R^{\infty * m}
function A = Amatrix(n , m , F) 
% F are the collocation of coefficient function 
% F = {chebfun(a^0(x)), chebfun(a^1(x)), ... ,chebfun(a^(m-1)(x))}
    A = zeros(n , m);
    x = chebfun(@(x)x);
    for i = 1 : 1 : m
        for k = 1 : 1 : i
            a = chebcoeffs(factorial(i-1) / factorial(i-k) * F{k} * x^(i-k));
            A(1 : length(a) , i) = A(1 : length(a) , i) + a;
        end
    end
end

function B = HighOrderBC(n , p , point)
% high order boundary condition of US
% n is length of chebyshev polynomials
% p is order of differential of chebyshev polynomials
% point = 1 or -1
    B = zeros(1 , n);
    for i = 1 : 1 : n
        b = (point)^(i + p - 1);
        for k = 0 : 1 : p - 1
           b = b * ((i-1)^2 - k^2) / (2 * k + 1);
        end
        B(1 , i) = b;
    end    
end

function x = upperTriangularSolve(U, b)
    % upperTriangularSolve 解决上三角线性方程组 Ux = b
    % U 是一个上三角矩阵
    % b 是一个向量
    % 返回解向量 x

    n = size(U, 1); % 获取矩阵的行数
    x = zeros(n, 1); % 初始化解向量

    for i = n:-1:1
        sum = b(i);
        for j = (i+1):n
            sum = sum - U(i,j) * x(j);
        end
        x(i) = sum / U(i,i);
    end
end


