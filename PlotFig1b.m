clear
clc
lw = 1;
fz = 20;

%% plot Fig1b form Julia
% read date
data = readmatrix('Fig1b.txt');
e = data(:, 1); % u64 - uexact
n = data(:, 2); % n

% plot 
semilogy(n , e, 'k', 'LineWidth',lw)

% legend('transport','heat', 'Location','east')
xlabel('$n$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\|\varepsilon\|_w$', 'Interpreter', 'latex', 'FontSize', 14)
set(gcf , 'Color' , 'w')
axis([0 350 1e-16 1e-0])
% set(gca,'yTick',[1e-21 1e-20 1e-19 1e-18 1e-17 1e-16 1e-15 1e-14]) 
% set(gca,'xTick',[0 10000 20000 30000 40000 50000]) 
% pbaspect([16 9 1])

exportgraphics(gcf, ['Float64-Uexact','.pdf'], 'Resolution', 300)