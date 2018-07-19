load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbiFinal.mat')
UbiSL = SL;
load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataIgorFinal.mat')
IgorSL = SL;
load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataKatoFinal.mat')
KatoSL = SL;

ContraRT = []; IpsiRT = [];
for S = 1:3
    switch S
        case 1
            SL = UbiSL;
        cs = [.5,.8,.5];
            
        case 2
            SL = IgorSL;
            cs = [.8,.5,.5];
            
        case 3
            SL = KatoSL;
            cs = [.5,.5,.8];
            
    end
    for i = 1:length(SL)
        if (~isempty(SL(i).Bad) || strcmp(SL(i).Condition, 'Control')...
                || strcmp(SL(i).Condition, 'NaN') || strcmp(SL(i).Condition, 'nostim')...
                || isempty(SL(i).trig1) || strcmp(SL(i).Condition(end),'R'))
            continue;
        end
        
        if(length(SL(i).Condition)<6 || ~strcmp(SL(i).Condition(1:6),'Contra') || SL(i).NormDelay<-100)
            continue;
        end
        
        if(strcmp(SL(i).StimHemi,'L'))
            CTrials = SL(i).righttrials;
            CRT = SL(i).rts_r;
            ITrials = SL(i).lefttrials;
            IRT = SL(i).rts_l;
        else
            ITrials = SL(i).righttrials;
            IRT = SL(i).rts_r;
            CTrials = SL(i).lefttrials;
            CRT = SL(i).rts_l;
        end
        
        bounds = [find(CTrials(:,2)<SL(i).trig1(1),1,'last'),find(CTrials(:,1)>SL(i).trig1(end),1)];
        
        ContraRT(end+1,1) = nanmedian(CRT(1:bounds(1)));
        ContraRT(end,2) = nanmedian(CRT(bounds(1)+1:bounds(2)-1));
        ContraRT(end,3) = nanmedian(CRT(bounds(2):end));
        
        bounds = [find(ITrials(:,2)<SL(i).trig1(1),1,'last'),find(ITrials(:,1)>SL(i).trig1(end),1)];
        
        IpsiRT(end+1,1) = nanmedian(IRT(1:bounds(1)));
        IpsiRT(end,2) = nanmedian(IRT(bounds(1)+1:bounds(2)-1));
        IpsiRT(end,3) = nanmedian(IRT(bounds(2):end));
        
    end
end


subplot(1,2,1); hold on;
% plot(ContraRT','Color',cs)
% for i = 1:size(ContraRT,1)
%     hold on; scatter([1,2,3],ContraRT(i,:),'MarkerEdgeColor',cs);
% end
b = bar([1,2,3],mean(ContraRT),1,'w');
uistack(b,'bottom');
set(gca,'XTick',1:3);
set(gca,'XTickLabel',{'PRE','STIM','POST'});
title('Contra Trials'); ylim([200,500]); xlim([0.5,3.5]);
ylabel('RT (ms)')
groups={[1,2],[2,3],[1,3]};  
[h,p] = ttest(ContraRT(:,2),ContraRT(:,1));
stats(1) = p;
[h,p] = ttest(ContraRT(:,2),ContraRT(:,3));
stats(2) = p;
[h,p] = ttest(ContraRT(:,1),ContraRT(:,3));
stats(3) = p;
u.sigstar(groups,stats);
ylim([200,550]);

subplot(1,2,2); hold on;
% plot(IpsiRT','Color',cs)
% for i = 1:size(IpsiRT,1)
%     hold on; scatter([1,2,3],IpsiRT(i,:),'MarkerEdgeColor',cs);
% end
b = bar([1,2,3],mean(IpsiRT),1,'w');
uistack(b,'bottom');
set(gca,'XTick',1:3);
set(gca,'XTickLabel',{'PRE','STIM','POST'});
title('Ipsi Trials'); xlim([0.5,3.5])
groups={[1,2],[2,3]};  stats = [];
[h,p] = ttest(IpsiRT(:,2),IpsiRT(:,1));
stats(1) = p;
[h,p] = ttest(IpsiRT(:,2),IpsiRT(:,3));
stats(2) = p;
u.sigstar(groups,stats);
ylim([200,550]);

I = plot(nan,nan,'color', [.8,.5,.5],'Linewidth',1.5); hold on;
K = plot(nan,nan,'color', [.5,.5,.8],'Linewidth',1.5); hold on;
U = plot(nan,nan,'color', [.5,.8,.5],'Linewidth',1.5);
legend([I,K,U],{'Monkey I','Monkey K','Monkey U'}); legend boxoff;


