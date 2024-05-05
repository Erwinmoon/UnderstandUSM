clear
clc
lw = 1;
fz = 20;
%% 
% test for 10-order function 

n = 10 : 10 : 500;
err = zeros(1 , length(n));

for i = 1 : 1 : length(n)
    n1 = n(i);
    n2 = ceil(1.01 * n1);

    u1 = Example10order(n1);
    u2 = Example10order(n2);
    
    e = u1 - u2(1 : n1);

    err(i) = norm(e , 2);
end

% plot 
semilogy(n , err, 'k', 'LineWidth',lw)

% legend('transport','heat', 'Location','east')
xlabel('$n$', 'Interpreter', 'latex', 'FontSize', 14)
% ylabel('Cauchy error', 'Interpreter', 'latex', 'FontSize', 14)
set(gcf , 'Color' , 'w')
axis([0 n(end) 1e-330 1e-0])
% exportgraphics(gcf, ['10orderCauchy','.pdf'], 'Resolution', 300)

function u = Example10order(n)
% example for 10-order example

    p = 4;
    f = zeros(n , 1);
    S0 = convertmat(n, 0 , 0);
    S1 = convertmat(n, 1 , 1);
    S2 = convertmat(n, 2 , 2);
    S3 = convertmat(n, 3 , 3);
    S4 = convertmat(n, 4 , 4);
    S5 = convertmat(n, 5 , 5);
    S6 = convertmat(n, 6 , 6);
    S7 = convertmat(n, 7 , 7);
    S8 = convertmat(n, 8 , 8);
    S9 = convertmat(n, 9 , 9);
    S10 = convertmat(n, 10 , 10);
    
    M8 = ultraS.multmat(n , chebfun(@(x)cosh(x)) , 8);
    M6 = ultraS.multmat(n , chebfun(@(x)x^2) , 6);
    M4 = ultraS.multmat(n , chebfun(@(x)x^4) , 4);
    M2 = ultraS.multmat(n , chebfun(@(x)cos(x)) , 2);
    M0 = ultraS.multmat(n , chebfun(@(x)x^2) , 2);
    
    D10 = diffmat(n , 10);
    D8 = diffmat(n , 8);
    D6 = diffmat(n , 6);
    D4 = diffmat(n , 4);
    D2 = diffmat(n , 2);
    D0 = diffmat(n , 0);
    

    l = -1;
    r = 1;
    B = zeros(2 * p , n);
    for i = 1 : 1 : p + 1
       B(2 * i - 1 , :) = HighOrderBC(n , i - 1 , l);
    end
    for i = 1 : 1 : p + 1
       B(2 * i , :) = HighOrderBC(n , i - 1 , r);
    end

    % LHS
    L = zeros(n , n);
    L = D10 + S9 * S8 *  M8 * D8 + S9 * S8 * S7 * S6 * M6 * D6 + ...
        S9 * S8 * S7 * S6 * S5 * S4 *  M4 * D4 + ...
        S9 * S8 * S7 * S6 * S5 * S4 *  S3 * S2 * M2 * D2 + ...
        S9 * S8 * S7 * S6 * S5 * S4 *  S3 * S2 * S1 * S0 * M0 * D0;
    L(2 * (p + 1) + 1 : end , :) = L(1 : n - 2 * (p + 1) , :);
    L(1 : 2 * (p + 1) , :) = B;

    % RHS
    f = S10 * f;
    f(2 * (p + 1) + 1 : end) = f(1 : n - 2 * (p + 1));
    f(1) = 0;
    f(2) = 0;
    f(3) = 1;
    f(4) = 1;
    f(5) = 0;
    f(6) = 0;
    f(7) = 0;
    f(8) = 0;
    f(9) = 0;
    f(10) = 0;

    % solve
    L = full(L);
    % [Q , R] = qr(L,0);
    [Q,R] = givensQR(L,30);
    u = upperTriangularSolve(R , Q' * f);
    % u = R(1:n,end);
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

function [Q, R] = givensQR(A,m)
    [n, ~] = size(A);
    Q = eye(n); 
    R = A;

    for j = 1:n
        for i = min(n, j + m):-1:(j + 1)
            if R(i, j) ~= 0
                [c, s] = givensCoefficients(R(i-1, j), R(i, j));
                
                % 更新R
                G = [c, s; -s, c];
                R([i-1, i], j:min(j+m, n)) = G' * R([i-1, i], j:min(j+m, n));
                
                % 更新Q
                Q(:, [i-1, i]) = Q(:, [i-1, i]) * G;
            end
        end
    end
end

function [c, s] = givensCoefficients(a, b)
    if b == 0
        c = 1;
        s = 0;
    else
        if abs(b) > abs(a)
            t = -a / b;
            s = 1 / sqrt(1 + t^2);
            c = s * t;
        else
            t = -b / a;
            c = 1 / sqrt(1 + t^2);
            s = c * t;
        end
    end
end