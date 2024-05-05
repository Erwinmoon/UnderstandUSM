clear
clc
lw = 1;
fz = 20;

%% plot Fig3 form Julia
% read date
data = readmatrix('Fig3.txt');
e = data(:, 1); % u64-ubf
mybound = data(:, 2); % bound by our condition number
cbound = data(: , 3); % bound by \kappa
n = data(:, 4); % n

% plot 
semilogy(n , e, 'k', 'LineWidth',lw)
hold on
semilogy(n , mybound, '--k', 'LineWidth',lw)
semilogy(n , cbound, '-.k', 'LineWidth',lw)
hold off


lgd = legend('$\varepsilon_S$','$\mathrm{cond}_{E,b}\epsilon_{mach}$','$\kappa_{2}\epsilon_{mach}$','Interpreter', 'latex', 'Location','east');
set(lgd, 'Position', [0.6, 0.5, 0.2, 0.2]);
set(lgd, 'FontSize', fz); 
xlabel('$n$', 'Interpreter', 'latex', 'FontSize', 14)
% ylabel('Cauchy error', 'Interpreter', 'latex')
set(gcf , 'Color' , 'w')
% axis([0 400 10^(-16.3) 10^(-15.8)])
% set(gca,'yTick',[1e-21 1e-20 1e-19 1e-18 1e-17 1e-16 1e-15 1e-14]) 
% set(gca,'xTick',[0 10000 20000 30000 40000 50000]) 
% pbaspect([16 9 1])

exportgraphics(gcf, ['eps1e-2','.pdf'], 'Resolution', 300)