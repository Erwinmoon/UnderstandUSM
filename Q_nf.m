clear
clc

hold off, offset = 0;
subplot('position',[.10 .25 .70 .70]);
LW = 'Linewidth'; lw = 0.7;
MS = 'MarkerSize'; ms = 12;
FS = 'FontSize'; fs = 15;
IN = 'interpret'; LX = 'latex';
L = 12; m = 6; m2 = m/2;
close all
set(gcf , 'Color' , 'w')
figure(1) 
hold on

% Qnf
text(offset+.3,-m2,'$\bf Q_nf$',FS,fs,IN,LX)

% Equals sign
offset = offset + 2.5;
plot(offset+[.2 .8],(-m2+.1)*[1 1],'k',LW,lw)
plot(offset+[.2 .8],(-m2-.1)*[1 1],'k',LW,lw)

% Q_nf vector
offset = offset + 1.5;
Qnf = fill(offset+[0 0 1 1],[1 6 6 1]-m,.69*[1 1 1]); set(Qnf, 'EdgeColor', 'None'),
plot(offset+[0 0 1 1 0],[1 6 6 1 1]-m,'k',LW,lw)
axis([0 30 -17 0]),axis square, axis off

% ,
offset = offset + 1.5;
text(offset+.1,-m2-.1,'$,$',FS,fs,IN,LX)

% Q_{n+s}f^{n+s}
offset = offset + 3;
text(offset+.3,-m2,'$\bf Q_{n+s}f^{n+s}$',FS,fs,IN,LX)

% Equals sign
offset = offset + 4.6;
plot(offset+[.2 .8],(-m2+.1)*[1 1],'k',LW,lw)
plot(offset+[.2 .8],(-m2-.1)*[1 1],'k',LW,lw)

% Q_{n+s}f^{n+s} vector
offset = offset + 1.5;
Qnf = fill(offset+[0 0 1 1],[0 6 6 0]-m,.69*[1 1 1]); set(Qnf, 'EdgeColor', 'None'),
plot(offset+[0 0 1 1 0],[0 6 6 0 0]-m,'k',LW,lw)

% Equals sign
offset = offset + 1.5;
plot(offset+[.2 .8],(-m2+.1)*[1 1],'k',LW,lw)
plot(offset+[.2 .8],(-m2-.1)*[1 1],'k',LW,lw)

% \phi_1
offset = offset + 1.5;
Qnf = fill(offset+[0 0 1 1],[2 6 6 2]-m,.69*[1 1 1]); set(Qnf, 'EdgeColor', 'None'),
plot(offset+[0 0 1 1 0],[0 6 6 0 0]-m,'k',LW,lw)

% RightPlus
offset = offset + 1 + 0.5;
plot(offset+[.2 .8],(-m2)*[1 1],'k',LW,lw)
plot(offset+[.5 .5],[-m2+0.2 -m2-0.2],'k',LW,lw)

% \phi_2
offset = offset + 1.5;
Qnf = fill(offset+[0 0 1 1],[0 2 2 0]-m,.69*[1 1 1]); set(Qnf, 'EdgeColor', 'None'),
plot(offset+[0 0 1 1 0],[0 6 6 0 0]-m,'k',LW,lw)

% exportgraphics(gcf, 'Qnf.pdf', 'Resolution', 300);
