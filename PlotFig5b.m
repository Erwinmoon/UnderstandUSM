clear
clc
lw = 1;
fz = 20;

%% plot Fig5b form Julia
% read date
data = readmatrix('Fig5b.txt');
c = data(:, 1); % cauchy
e = data(:, 2); % u64 - uexact
n = data(:, 3); % n

% plot 
semilogy(n , c, 'k', 'LineWidth',lw)
hold on
semilogy(n , e, '--k', 'LineWidth',lw)
hold off


lgd = legend('$Cauchy \; error$','$\|\varepsilon\|_w$','Interpreter', 'latex')
xlabel('$n$', 'Interpreter', 'latex', 'FontSize', 14)
set(lgd, 'FontSize', fz); 
% ylabel('Cauchy error', 'Interpreter', 'latex')
set(gcf , 'Color' , 'w')
% axis([0 400 10^(-16.3) 10^(-15.8)])
% set(gca,'yTick',[1e-21 1e-20 1e-19 1e-18 1e-17 1e-16 1e-15 1e-14]) 
% set(gca,'xTick',[0 10000 20000 30000 40000 50000]) 
% pbaspect([16 9 1])

exportgraphics(gcf, ['1e-5','.pdf'], 'Resolution', 300)