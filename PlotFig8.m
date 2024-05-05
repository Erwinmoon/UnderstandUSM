clear
clc
lw = 1;
fz = 14;

%% plot Fig8 form Julia
% read date
data = readmatrix('Fig8.txt');
e1 = data(:, 1); % u64 - uexact 1e-2
e2 = data(:, 2); % u64 - uexact 1e-5
n = data(:, 3); % n

% plot 
semilogy(n , e1, 'k', 'LineWidth',lw)
hold on
semilogy(n , e2, '-.k', 'LineWidth',lw)
hold off

text(250, 1e-14, '$\mu = 10^{-2}$', 'Interpreter', 'latex', 'FontSize', 12);
text(560, 10^(-10.7), '$\mu = 10^{-5}$', 'Interpreter', 'latex', 'FontSize', 12);

% lgd = legend('$\|\varepsilon\|_w(1e-2)$','$\|\varepsilon\|_w(1e-5)$','Interpreter', 'latex');
xlabel('$n$', 'Interpreter', 'latex', 'FontSize', 14)
% set(lgd, 'FontSize', fz); 
% set(lgd, 'Location', 'southwest'); 
ylabel('$\|\varepsilon\|_w$', 'Interpreter', 'latex')
set(gcf , 'Color' , 'w')
axis([0 1000 10^(-16) 10^(2)])
% set(gca,'yTick',[1e-21 1e-20 1e-19 1e-18 1e-17 1e-16 1e-15 1e-14]) 
% set(gca,'xTick',[0 10000 20000 30000 40000 50000]) 
% pbaspect([16 9 1])

exportgraphics(gcf, ['veps_1e-2_1e-5','.pdf'], 'Resolution', 300)