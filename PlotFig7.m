clear
clc
lw = 1;
fz = 9;

%% plot Fig7 form Julia
% read date
data = readmatrix('Fig7.txt');
c1 = data(:, 1); % cauchy for 1e-2
Q1 = data(:, 2); % Q(1:2,end-2,end) 1e-2
c2 = data(:, 3); % cauchy for 1e-5
Q2 = data(:, 4); % Q(1:2,end-2,end) 1e-5
n = data(:, 5); % n

k1 = 0;
k2 = 0;

for i = 1 : 1 : length(Q1)
    if abs(c1(i)) +  abs(Q1(i)) ~= 0
        k1 = i;
    end
    if abs(c2(i)) +  abs(Q2(i)) ~= 0
        k2 = i;
    end
end

% plot 
figure(1)
hold on
plot(n(1 : k1) , log10(c1(1 : k1)), 'k', 'LineWidth',lw)
plot(n(1 : k1) , log10(Q1(1 : k1)), '--k', 'LineWidth',lw)
plot(n(1 : k2) , log10(c2(1 : k2)), 'k', 'LineWidth',lw)
plot(n(1 : k2) , log10(Q2(1 : k2)), '--k', 'LineWidth',lw)
hold off

text(250, -200, '$\mu = 10^{-2}$', 'Interpreter', 'latex', 'FontSize', 12);
text(560, -100, '$\mu = 10^{-5}$', 'Interpreter', 'latex', 'FontSize', 12);

lgd = legend('$Cauchy \; error$','$\|Q_n$\texttt(1,n-3:n)$\|_2 + \|Q_n$\texttt(2,n-3:n)$\|_2$','Interpreter', 'latex');
xlabel('$n$', 'Interpreter', 'latex', 'FontSize', 14)
set(lgd, 'FontSize', fz); 
set(lgd, 'Location', 'northeast'); 
% ylabel('Cauchy error', 'Interpreter', 'latex')
set(gcf , 'Color' , 'w')
axis([0 1000 -350 10^(1.5)])
% 自定义 x 轴刻度
yticks([-350 -300 -250 -200 -150 -100 -50 0]); % 设置 x 轴的刻度值
% 自定义 x 轴刻度的标签文本
yticklabels({'$10^{-350}$', '$10^{-300}$', '$10^{-250}$', '$10^{-200}$', '$10^{-150}$', '$10^{-100}$', '$10^{-50}$', '$10^{0}$'});
set(gca, 'TickLabelInterpreter', 'latex')
% set(gca,'yTick',[1e-350 1e-300 1e-200 1e-100 1e-0]) 
% set(gca,'xTick',[0 10000 20000 30000 40000 50000]) 
% pbaspect([16 9 1])

exportgraphics(gcf, ['CauchyQ','.pdf'], 'Resolution', 300)