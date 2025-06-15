clear;close all;clc;FS=12;
%% List of file names
lst=dir('.\micRawData\');

lst=lst(3:end); % Exclude markers of directories (dots)
for N=1:length(lst);
    space1=strfind(lst(N).name,' ');
    nameN0=lst(N).name(1:space1(1)-1); 
    nameN(N)=sum(double(nameN0));
    nameN=nameN';
end
%%
for N=1:length(lst);
    space1=strfind(lst(N).name,' ');
    nameN0=lst(N).name(1:space1(1)-1); 
    nameN(N)=sum(double(nameN0));
    nameN=nameN';
end
%%
% Read data from a file
k=1;
clr={'blue','red','magenta','green'};
for N=7:9;%1:length(lst);
    Tbl=xlsread(['.\micRawData\',lst(N).name]);
    space1=strfind(lst(N).name,' ');
    name=lst(N).name(1:space1(1)-1);
    % Fluorescence data
    % Actual concentrations
    conc=Tbl(10,1:9);
    % Fluorescence witout the drug
    data0=Tbl(2:9,10)';
    data0r=repmat(data0,8,1);
    data0rs=reshape(data0r,64,1);
    data0rsr=repmat(data0rs,1,9);
    % Fluorescence data with the drug
    data=Tbl(2:9,1:9);
    datar=repmat(data,8,1);
    % Normed data
    F=datar./data0rsr;
    Q5=quantile(F,0.05);
    Q95=quantile(F,0.95);
%     plot(conc,Q,'o:');
    indB=find((Q5>0.1)&(Q95>0.1));indB=indB(1);
    indT=find((Q5<0.1)&(Q95<0.1));indT=indT(end);
    FPint=[];cint=[];
    for j=indT:indB
        ind=(F(:,j)>0)&(F(:,j)<1);
        FP=log((1-F(ind,j))./F(ind,j))/log(9)-1;
%         FPQ5(j)=quantile(FP,0.05);
%         FPQ95(j)=quantile(FP,0.95);
        FPint=[FPint;FP];
        cint=[cint;repmat(conc(j),length(FP),1)];
    end
    % 
      plot(conc,log((1-F)./F)/log(9)-1,'o','color',clr{k},'MarkerSize',4)
      hold on
      plot(cint,FPint,'x','color',clr{k},'MarkerSize',8)
      hold on
    
    set(gca,'XScale','log')
    set(gca,'XTick',[0.01 0.02 0.04 0.08 0.16 0.32 0.64 1.28 2.56 5.12 10.24 20.48])
    box on
    xlim([0.005 25])    
    hold on
    plot([0.005 40],[0 0],':','color','black')
    logcmod=fit(FPint,log(cint),'poly1','Robust','bisquare');
    MICint(N,:)=exp(predint(logcmod,0));
    pint=predint(logcmod,[-1:0.1:1]);
    plot(exp(logcmod(FPint)),FPint,'color',clr{k},'LineWidth',1.5)
    plot(exp(pint),[-1:0.1:1],'--','color',clr{k},'LineWidth',1.5)    
    k=k+1;
end
    xlim([0.04 40])
    xlabel('Concentration, µg/ml')
    ylabel('Scaled Fisher-Pry plot')
%    title(name)
    set(gca,'FontSize',FS,'FontName','Times')
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 18 9],'PaperSize',[18 9])
print('FigExample2','-dpng','-r300','-f')