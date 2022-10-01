Nt = size(Timepoint,2);

figure('visible','off')
y = [1:Nt];
x = -.5*ones(size(y));
scatter(x,y,60,C,'s','filled')
hold on
for ii = y
    text(.5,ii,Timepoint{ii},'FontSize',default_fs)
end
y = Nt+3:-1:Nt+1;
for ct=1:3
    scatter(x(ct),y(ct),50,'k',Markers{ct})
    text(.5,y(ct),CellTypes{ct},'FontSize',default_fs)
end
axis off

dim = [2 6];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto','PaperPosition',[0  0 dim],'PaperSize',[dim]);
print(gcf,'Fig/legend','-dpdf');


