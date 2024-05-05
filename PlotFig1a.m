clear
clc
lw = 1;
fz = 20;

%% plot Fig1a form Julia
% read date
data = readmatrix('Fig1a.txt');
c = data(:, 1); % conditon number
n = data(:, 2); % n

% plot 
plot(n , c, 'k', 'LineWidth',lw)

% legend('transport','heat', 'Location','east')
xlabel('$n$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\kappa_{\infty}$', 'Interpreter', 'latex', 'FontSize', 14)
set(gcf , 'Color' , 'w')
% axis([0 50000 1e-21 1e-14])
% set(gca,'yTick',[1e-21 1e-20 1e-19 1e-18 1e-17 1e-16 1e-15 1e-14]) 
% set(gca,'xTick',[0 10000 20000 30000 40000 50000]) 
% pbaspect([16 9 1])

exportgraphics(gcf, ['Condition','.pdf'], 'Resolution', 300)