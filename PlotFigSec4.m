clear
clc
lw = 1;
fz = 14;

%% plot Fig1a form Julia
% read date
data1 = readmatrix('shenfun.txt');
data2 = readmatrix('shenfun2.txt');
data3 = readmatrix('kuidu.txt');
n = data1(:, 1); % n
errs = data2; % shenfun error
errd = data3(: , 2); % kuidu error
errcauchys = data1(:, 2);
errcauchyd = data3(: , 3);


% plot 
figure(1)
semilogy(n , errs, 'k', 'LineWidth',lw)
hold on
semilogy(n , errd, '--k', 'LineWidth',lw)
hold off

lgd = legend('$BG$','$IR \; tau$','Interpreter', 'latex', 'Location','east');
set(lgd, 'Position', [0.65, 0.5, 0.2, 0.2]);
set(lgd, 'FontSize', fz); 
xlabel('$n$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\|\varepsilon\|$', 'Interpreter', 'latex', 'FontSize', 14)
axis([0 350 1e-16 1e-0])
set(gcf , 'Color' , 'w')
exportgraphics(gcf, 'sec4error.pdf', 'Resolution', 300);

figure(2)
semilogy(n , errcauchys, 'k', 'LineWidth',lw)
hold on
semilogy(n , errcauchyd, '--k', 'LineWidth',lw)
hold off

lgd = legend('$BG$','$IR \; tau$','Interpreter', 'latex', 'Location','east');
set(lgd, 'Position', [0.65, 0.5, 0.2, 0.2]);
set(lgd, 'FontSize', fz); 
xlabel('$n$', 'Interpreter', 'latex', 'FontSize', 14)
% ylabel('$\varepsilon_C$', 'Interpreter', 'latex', 'FontSize', 14)
set(gcf , 'Color' , 'w')
exportgraphics(gcf, 'sec4cauchyerror.pdf', 'Resolution', 300);

% axis([0 50000 1e-21 1e-14])
% set(gca,'yTick',[1e-21 1e-20 1e-19 1e-18 1e-17 1e-16 1e-15 1e-14]) 
% set(gca,'xTick',[0 10000 20000 30000 40000 50000]) 
% pbaspect([16 9 1])

