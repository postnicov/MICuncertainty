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
for N=1:length(lst);
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
        FPQ5(j)=quantile(FP,0.05);
        FPQ95(j)=quantile(FP,0.95);
        FPint=[FPint;FP];
        cint=[cint;repmat(conc(j),length(FP),1)];
    end
    logcmod=fit(FPint,log(cint),'poly1','Robust','bisquare');
    MICint(N,:)=exp(predint(logcmod,0));
end

figure
k=1;
kk=1;
Nd=length(lst);
clr={'blue','red','magenta','green'};
mrk={'o','s','d','x'};
mrk={'.','.','.','.'};
for N=1:Nd;
    plot(MICint(N,:),[k+kk/20 k+kk/20],'color',clr{kk},'marker',mrk{kk},'LineWidth',2);
    hold on
    kk=kk+1;
    space1=strfind(lst(N).name,' ');
    Yname{k}=lst(N).name(1:space1(1)-1);
    if (N<Nd)&(nameN(N+1)~=nameN(N))
        k=k+1;
        kk=1;
    end
end
%
load TblLit
for j=1:10
    if sum(TblLit(j,1:2))>0
        plot(TblLit(j,1:2),[j-1/20 j-1/20],'-','color','black','LineWidth',1.5);
    end
    if sum(TblLit(j,3:4))>0
        plot(TblLit(j,3:4),[j-2/20 j-2/20],'--','color','black','LineWidth',1.5);
    end
    if sum(TblLit(j,5:6))>0
        plot(TblLit(j,5:6),[j-3/20 j-3/20],'-.','color','black','LineWidth',1.5);
    end    
end
%
    set(gca,'XScale','log')
    set(gca,'FontSize',FS,'FontName','Times')
    set(gca,'XTick',[0.01 0.02 0.04 0.08 0.16 0.32 0.64 1.28 2.56 5.12 10.24 20.48])
    box on
    xlim([0.04 10.24]) 
    LY=length(Yname);
    xlabel('Concentration, µg/ml')
    set(gca,'YTick',[1:LY],'YTicklabel',Yname);
    ylim([0.5 LY+0.5])
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 15 18],'PaperSize',[15 18])

print('FigRanges','-dpng','-r300')