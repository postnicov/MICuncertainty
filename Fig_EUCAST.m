clear;close all;clc;FS=11;

Tbl=readtable('EUCAST_MIC_30_05_2025.xlsx');

data=table2array(Tbl(2:end,2:end));
data(9,9)=0;
[N,L]=size(data);
conc=table2array(Tbl(1,2:end));

for j=1:N;
    subplot(N,1,j)
    S=sum(data(j,:));
    P=100*data(j,:)/S;
    bar(P)
    ylim([0 60])
    set(gca,'YTick',[0 30 60],'YTicklabel',['0%',Tbl.Var1(j+1),'60%'])
    set(gca,'XTicklabel',conc)
    set(gca,'FontSize',FS,'FontName','Times');
end
xlabel('Concentration, µg/ml')
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 15 18],'PaperSize',[15 18])
print('FigEUCAST','-dpng','-r150')