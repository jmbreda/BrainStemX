function Phases_polar_plot()
%clear all; close all; clc;

set(0,'DefaultTextFontSize',50)
% Saturation = r
% Hue = phi/2*pi

R = linspace(.8,1,8);
PHI = linspace(0,2*pi,200);
PHI = PHI(1:end-1);

x = [];
y = [];
col = [];
S = [];
for r = R
    for phi = PHI
        x = [x; r*sin(phi)];
        y = [y; r*cos(phi)];
        col = [col; hsv2rgb([phi/(2*pi) r 1])];
        S = [S; 20*r];
    end
end
disp('plot')
figure('visible','off');
scatter(x,y,S.^2,col,'filled')
hold on;
r = .75;
Phi = [0 2*pi/3 4*pi/3];
Phase = {'NSC','BP','NBN'};
k = 1;
text(r*sin(Phi(k)),r*cos(Phi(k)),Phase{k},'HorizontalAlignment','center','VerticalAlignment','top');
k = 2;
text(r*sin(Phi(k)),r*cos(Phi(k)),Phase{k},'HorizontalAlignment','right','VerticalAlignment','bottom');
k = 3;
text(r*sin(Phi(k)),r*cos(Phi(k)),Phase{k},'HorizontalAlignment','left','VerticalAlignment','bottom');
axis off
dim = [32 32];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto','PaperPosition',[0 0 dim],'PaperSize',[dim]);
print(gcf,'Fig/Celltype_phase_plot_s','-dpdf');
print(gcf,'Fig/Celltype_phase_plot_s','-dpng');

figure('visible','off');
scatter(x,y,S.^2,col,'filled')
hold on;
r = .75;
Phi = [0 2*pi/3 4*pi/3];
Phase = {'E','N','G'};
k = 1;
text(r*sin(Phi(k)),r*cos(Phi(k)),Phase{k},'HorizontalAlignment','center','VerticalAlignment','top');
k = 2;
text(r*sin(Phi(k)),r*cos(Phi(k)),Phase{k},'HorizontalAlignment','right','VerticalAlignment','bottom');
k = 3;
text(r*sin(Phi(k)),r*cos(Phi(k)),Phase{k},'HorizontalAlignment','left','VerticalAlignment','bottom');
axis off
dim = [32 32];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto','PaperPosition',[0 0 dim],'PaperSize',[dim]);
print(gcf,'Fig/NSC_phase_plot_s','-dpdf');
print(gcf,'Fig/NSC_phase_plot_s','-dpng');


