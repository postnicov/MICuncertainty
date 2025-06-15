clear;close all;clc;FS=12;
%% List of file names
lst=dir('.\micRawData\');
fid = fopen('TableOut.txt', 'wt');
fprintf(fid, '%s\n', '\begin{longtable}{cccccccccc}');
fprintf(fid, '%s\n', '\hline');
fprintf(fid, '%s\n', 'No&Name&$f_{min}$&$CI_{f_{min}}$&$f_{max}$&$CI_{f_{max}}$&$IC50$&$CI_{IC50}$&$\alpha$&$CI_{\alpha}$\\');
fprintf(fid, '%s\n', '\hline');
lst=lst(3:end); % Exclude markers of directories (dots)
for N=1:length(lst);
    space1=strfind(lst(N).name,' ');
    nameN0=lst(N).name(1:space1(1)-1); 
    nameN(N)=sum(double(nameN0));
    nameN=nameN';
end
% Read data from a file
for N=1:length(lst);
    Tbl=xlsread(['.\micRawData\',lst(N).name]);
    space1=strfind(lst(N).name,' ');
    name=lst(N).name(1:space1(1)-1);
    % Fluorescence data
    data=Tbl(2:9,1:10);
    data=data/median(data(:,10));
    conc=Tbl(10,1:9);
    Mdata=median(data(:,1:9));
    sigma=1.4826*mad(data(:,1:9),1);
    %%
    errorbar(conc,Mdata,sigma,'o','color','black');
    hold on
    %%
    Hill=@(mn,mx,IC50,alpha,x)  mn+(mx-mn)./(1+(x/IC50).^alpha);
    b0=[min(Mdata) 1 mean(conc) 3];
    indL=find((Mdata<1.05));
    errorbar(conc(indL),Mdata(indL),sigma(indL),'o','color','blue');
    if length(indL)>3
    fitfun=fit(conc(indL)',Mdata(indL)',Hill,...
        'Lower',[0 -Inf 0 0],...
        'StartPoint', b0,'Robust','LAR');
    %%
    concx=logspace(-2,2,151);
    plot(concx,fitfun(concx),'color','red')
        set(gca,'XScale','log')
        set(gca,'XTick',[0.01 0.02 0.04 0.08 0.16 0.32 0.64 1.28 2.56 5.12 10.24 20.48])
        box on
        xlim([0.005 40])

    %
    [~,indmn]=min(Mdata);
    [~,indmx]=max(Mdata);
    %
    ylim([0 1.3])
    fname=['./plots/Fig',name,'.png'];
    %
    xlabel('Concentration, \mug/ml')
    ylabel('Normed fluorescence')
    title(name)
    set(gca,'FontSize',FS,'FontName','Times')
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 18 6],'PaperSize',[18 6])
    fname=['./plots/',num2str(N),'_',name,'.png'];
    print(fname,'-dpng','-r150','-f')

    ci=confint(fitfun);
    Sres=[num2str(N),'&',name,'&',...
          num2str(fitfun.mn,2),'&(',num2str(ci(1,1),3),',',num2str(ci(2,1),3),')&'...
          num2str(fitfun.mx,2),'&(',num2str(ci(1,2),3),',',num2str(ci(2,2),3),')&'...
          num2str(fitfun.IC50,2),'&(',num2str(ci(1,3),3),',',num2str(ci(2,3),3),')&'...
          num2str(fitfun.alpha,2),'&(',num2str(ci(1,4),3),',',num2str(ci(2,4),3),')\\'];
     fprintf(fid, '%s\n', Sres);
    end
    if (N<(length(lst))-1)&(nameN(N+1)~=nameN(N))
    close
    end
end
fprintf(fid, '%s\n', '\hline');
fprintf(fid, '%s\n', '\end{longtable}');
fclose(fid);

