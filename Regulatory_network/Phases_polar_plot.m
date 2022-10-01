function Phases_polar_plot()
%clear all; close all; clc;

set(0,'DefaultTextFontSize',50)
% Saturation = r
% Hue = phi/2*pi

R = linspace(.8,1,9);
%R = linspace(0.025,1,40);

x = [];
y = [];
col = [];
S = [];
for r = R
    N = 50+150*r;
    PHI = linspace(0,2*pi,N);
    PHI = PHI(1:end-1);
    for phi = PHI
        x = [x; r*sin(phi)];
        y = [y; r*cos(phi)];
        col = [col; hsv2rgb([phi/(2*pi) r 1])];
        S = [S; 100+100*r];
    end
end
disp('plot')
figure('visible','off');
scatter(x,y,S,col,'filled')

if 0
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
end
axis off
dim = [32 32];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto','PaperPosition',[0 0 dim],'PaperSize',[dim]);
print(gcf,'Fig/Phase_plot','-dpdf');
print(gcf,'Fig/Phase_plot','-dsvg');
close all

