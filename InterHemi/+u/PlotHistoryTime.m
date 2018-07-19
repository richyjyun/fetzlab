% load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbi.mat')
% SL = a.RemoveTwitch(SL); SL = u.AppendThreshold(SL); SL = a.AppendRT(SL); SL = a.AppendNormalizedDelay(SL); UbiSL = SL;
% load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataIgor.mat')
% SL = a.RemoveTwitch(SL); SL = u.AppendThreshold(SL); SL = a.AppendRT(SL); SL = a.AppendNormalizedDelay(SL); IgorSL = SL;
% load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataKato.mat')
% SL = a.RemoveTwitch(SL); SL = u.AppendThreshold(SL); SL = a.AppendRT(SL); SL = a.AppendNormalizedDelay(SL); KatoSL = SL;


ContraRT = {};
IpsiRT = {};
ContraTime = {};
IpsiTime = {};

PreContraRT = {};
PreContraTime = {};
PreIpsiRT = {};
PreIpsiTime = {};

PreContra = [];
PreIpsi = [];

for S = 1
    switch S
        case 1
            prevInd = 0;
            SL = UbiSL;
        case 2
%             prevInd = length(UbiSL);
            SL = IgorSL;
        case 3
            prevInd = length(UbiSL)+length(IgorSL);
            SL = KatoSL;
    end
    for i = 1:length(SL)
        if (~isempty(SL(i).Bad) || strcmp(SL(i).Condition, 'Control')...
                || strcmp(SL(i).Condition, 'NaN') || strcmp(SL(i).Condition, 'nostim')...
                || isempty(SL(i).trig1) || strcmp(SL(i).Condition(end),'R'))
            continue;
        end
        
        D = str2num(SL(i).Date);
        Bad = 20120515;%[20120508,20120515,20120612,20120814];
        if( any(Bad==D))
            continue;
        end

%         if(~(strcmp(SL(i).Condition(end),'M') || str2num(SL(i).Stim_Delay) > 250))
%             continue;
%         end
        if(~strcmp(SL(i).Condition(1:6),'Contra'))
            continue;
        end

        disp(SL(i).Date)
        
        %determine contra and ipsi trials
        Contra = []; Ipsi = []; CRT = []; IRT = [];
        if(strcmp(SL(i).StimHemi,'L'))
            Contra = SL(i).righttrials;
            Ipsi = SL(i).lefttrials;
            CRT = SL(i).rts_r;
            IRT = SL(i).rts_l;
        elseif(strcmp(SL(i).StimHemi,'R'))
            Contra = SL(i).lefttrials;
            Ipsi = SL(i).righttrials;
            CRT = SL(i).rts_l;
            IRT = SL(i).rts_r;
        end
        
        %determine which trials have stimulation
        StimTrials = Contra;
        Trials = [Contra(:,1);Ipsi(:,1)];
        RT = [CRT;IRT];
        Label = char([ones(length(CRT),1)*67;ones(length(IRT),1)*73]);  %uint16 code for C and I

        StimInd = zeros(length(StimTrials),1);
        for j = 1:length(SL(i).trig1)
            norm = StimTrials(:,1)-50-SL(i).trig1(j);
            ind = find(norm>0,1)-1;
            if(~isempty(ind))
                StimInd(ind) = 1;
            end
        end
        StimInd = [StimInd;zeros(length(Trials)-length(StimInd),1)];
        
        
        maxHist = 4;
        %set up order of all trials
        [Trials,Order] = sort(Trials);
        StimInd = StimInd(Order);
        Trim = [find(StimInd,1),find(StimInd,1,'last')+maxHist];
        StimInd = StimInd(Trim(1):Trim(2));
        PreTrials = Trials(1:Trim(1)-1,:);
        Trials = Trials(Trim(1):Trim(2));
        PreRT = RT;
        RT = RT(Order); PreRT = RT(1:Trim(1)-1); RT = RT(Trim(1):Trim(2));
        Label = Label(Order); PreLabel = Label(1:Trim(1)-1); Label = Label(Trim(1):Trim(2));
        ContraPre = nanmedian(PreRT(strfind(PreLabel','C')));
        IpsiPre = nanmedian(PreRT(strfind(PreLabel','I')));
        PreContra(end+1) = nanmedian(RT(strfind(Label','C')))-ContraPre;
        PreIpsi(end+1) = nanmedian(RT(strfind(Label','I')))-IpsiPre;
        
        %Store Pre history and counts. History being since last contra
        %trial.
        StimHist = nan(length(PreLabel),1); all = strfind(PreLabel','C');
        for j = 2:length(StimHist)
            last = all(find(all<j,1,'last'));
            if(~isempty(last))
                StimHist(j) = PreTrials(j,1)-PreTrials(last,1);
            end
        end
        inds = strfind(PreLabel','C');
        PreContraTime{end+1} = StimHist(inds);
        PreContraRT{end+1} = PreRT(inds)-ContraPre;
        
        inds = strfind(PreLabel','I');
        PreIpsiTime{end+1} = StimHist(inds);
        PreIpsiRT{end+1} = PreRT(inds)-IpsiPre;
            
        %Store Pre history and counts. History being since last contra
        %trial.
        StimHist = nan(length(Label),1);
        for j = 2:length(StimHist)
            last = find(StimInd(1:j-1) == 1,1,'last');
            if(~isempty(last))
                StimHist(j) = Trials(j,1)-Trials(last,1);
            end
        end
        
        inds = strfind(Label','C');
        ContraTime{end+1} = StimHist(inds);
        ContraRT{end+1} = RT(inds)-ContraPre;
        
        inds = strfind(Label','I');
        IpsiTime{end+1} = StimHist(inds);
        IpsiRT{end+1} = RT(inds)-IpsiPre;
        
    end
end

figure; subplot(2,2,1); maxTime = 15000; binWidth = 1000; edges = 0:binWidth:maxTime;
ContraAll = nan(length(ContraRT),length(edges)-1);
RT = []; Time = []; %fit = struct([]);
for(i = 1:length(ContraRT))
    hold on; 
%     [N,edges] = histcounts(ContraTime{i}(ContraTime{i}<maxTime));
%     bins = discretize(ContraTime{i},edges);
%     for j = 1:max(bins)
%         errorbar(edges(j)+binWidth/2,nanmean(ContraRT{i}(bins==j)),nanstd(ContraRT{i}(bins==i)/sqrt(sum(bins==j))),'ko')
%         ContraAll(i,j) = nanmean(ContraRT{i}(bins==j));
%     end
%     scatter(ContraTime{i},ContraRT{i},'k.');
    Time = [Time;ContraTime{i}(ContraTime{i}<maxTime)];
    RT = [RT;ContraRT{i}(ContraTime{i}<maxTime)];
%     tbl = table(ContraRT{i}(ContraTime{i}<maxTime),ContraTime{i}(ContraTime{i}<maxTime),'VariableNames',{'RT','Time'});
%     fit(end+1).f = fitlm(tbl,'RT~Time');
end
hold on; errorbar(edges(1:end-1)+binWidth/2,nanmean(ContraAll),nanstd(ContraAll)./sqrt(sum(~isnan(ContraAll))),'ro','MarkerSize',10);
% xlim([0,maxTime]); title('Contra Trials'); xlabel('Time since last stim (ms)'); ylabel('\Delta RT (ms)');
tbl = table(RT,Time/1000,'VariableNames',{'RT','Time'});
fit1 = fitlm(tbl,'RT~Time','RobustOpts','on');
h = plot(fit2); h(1).Marker = 'o'; h(1).Color = [0,0,0]; legend HIDE; h(2).LineWidth = 2; h(3).LineWidth = 1.5; h(4).LineWidth = 1.5;
title(['Contra Trials']);%. Slope: ',num2str(fit1.Coefficients{2,1}),' ms/s'])%,'; P-value: ',num2str(fit1.Coefficients{2,2})]); 
xlabel('Time since last stim (s)'); ylabel('\Delta RT (ms)','Interpreter','tex');
h = plot(nan,nan,'w'); legend(h,['Slope: ',num2str(fit1.Coefficients{2,1}),' ms/s']); legend('boxoff')
% hold on; xl = xlim; plot(xlim,[0,0],'k--')
% Intercept = []; Slope = [];
% for i = 1:length(fit)
%     Intercept(i,:) = [fit(i).f.Coefficients{1,1},fit(i).f.Coefficients{1,4}];
%     Slope(i,:) = [fit(i).f.Coefficients{2,1},fit(i).f.Coefficients{2,4}];
% end
% f1 = fit(Time(~isnan(RT)),RT(~isnan(RT)),'exp1'); hold on; plot(f);


RT = []; Time = [];
subplot(2,2,2); IpsiAll = nan(length(IpsiRT),length(edges)-1);
for(i = 1:length(IpsiRT))
    hold on; 
%     [N,edges] = histcounts(IpsiTime{i}(IpsiTime{i}<maxTime));
%     bins = discretize(IpsiTime{i},edges);
%     for j = 1:max(bins)
%         errorbar(edges(j)+binWidth/2,nanmean(IpsiRT{i}(bins==j)),nanstd(IpsiRT{i}(bins==i)/sqrt(sum(bins==j))),'ko')
%         IpsiAll(i,j) = nanmean(IpsiRT{i}(bins==j));
%     end
%     scatter(IpsiTime{i},IpsiRT{i},'k.');
    Time = [Time;IpsiTime{i}(IpsiTime{i}<maxTime)];
    RT = [RT;IpsiRT{i}(IpsiTime{i}<maxTime)];
end
hold on; errorbar(edges(1:end-1)+binWidth/2,nanmean(IpsiAll),nanstd(IpsiAll)./sqrt(sum(~isnan(IpsiAll))),'ro','MarkerSize',10);
% xlim([0,maxTime]); title('Ipsi Trials'); xlabel('Time since last stim (ms)'); ylabel('\Delta RT (ms)');
tbl = table(RT,Time/1000,'VariableNames',{'RT','Time'});
fit2 = fitlm(tbl,'RT~Time','RobustOpts','on');
h = plot(fit2); h(1).Marker = '.'; h(1).Color = [0,0,0]; legend HIDE; h(2).LineWidth = 2; h(4).LineWidth = 1.5;
title(['Ipsi Trials']);%. Slope: ',num2str(fit2.Coefficients{2,1}),' ms/s'])%,'; P-value: ',num2str(fit2.Coefficients{2,2})]); 
xlabel('Time since last stim (s)'); ylabel('\Delta RT (ms)','Interpreter','tex');
h = plot(nan,nan,'w'); legend(h,['Slope: ',num2str(fit1.Coefficients{2,1}),' ms/s']); legend('boxoff')
% hold on; xl = xlim; plot(xlim,[0,0],'k--')





% 
% RT = []; Time = []; %fit = struct([]);
% figure; PreContraAll = nan(length(PreContraRT),length(edges)-1);
% for(i = 1:length(PreContraRT))
%     hold on; 
% %     [N,edges] = histcounts(ContraTime{i}(ContraTime{i}<maxTime));
%     bins = discretize(PreContraTime{i},edges);
%     for j = 1:max(bins)
%         errorbar(edges(j)+binWidth/2,nanmean(PreContraRT{i}(bins==j)),nanstd(PreContraRT{i}(bins==i)/sqrt(sum(bins==j))),'ko')
%         PreContraAll(i,j) = nanmean(PreContraRT{i}(bins==j));
%     end
% %     scatter(PreContraTime{i},PreContraRT{i});
% %     Time = [Time;PreContraTime{i}(PreContraTime{i}<maxTime)];
% %     RT = [RT;PreContraRT{i}(PreContraTime{i}<maxTime)];
% %     tbl = table(ContraRT{i}(ContraTime{i}<maxTime),ContraTime{i}(ContraTime{i}<maxTime),'VariableNames',{'RT','Time'});
% %     fit(end+1).f = fitlm(tbl,'RT~Time');
% end
% hold on; errorbar(edges(1:end-1)+binWidth/2,nanmean(PreContraAll),nanstd(PreContraAll)./sqrt(sum(~isnan(PreContraAll))),'ro','MarkerSize',10);
% xlim([0,maxTime]); title('Pre Contra Trials'); xlabel('Time since last stim (ms)'); ylabel('\Delta RT (ms)');
% tbl = table(RT,Time,'VariableNames',{'RT','Time'});
% fitPre = fitlm(tbl,'RT~Time');
% 
% 
% RT = []; Time = []; %fit = struct([]);
% figure;
% for(i = 1:length(PreIpsiRT))
%     hold on; 
% %     [N,edges] = histcounts(ContraTime{i}(ContraTime{i}<maxTime));
% %     bins = discretize(ContraTime{i},edges);
% %     for j = 1:max(bins)
% %         errorbar(edges(j)+binWidth/2,nanmean(ContraRT{i}(bins==j)),nanstd(ContraRT{i}(bins==i)/sqrt(sum(bins==j))),'ko')
% %         ContraAll(i,j) = nanmean(ContraRT{i}(bins==j));
% %     end
%     scatter(PreIpsiTime{i},PreIpsiRT{i});
%     Time = [Time;PreIpsiTime{i}(PreIpsiTime{i}<maxTime)];
%     RT = [RT;PreIpsiRT{i}(PreIpsiTime{i}<maxTime)];
% %     tbl = table(ContraRT{i}(ContraTime{i}<maxTime),ContraTime{i}(ContraTime{i}<maxTime),'VariableNames',{'RT','Time'});
% %     fit(end+1).f = fitlm(tbl,'RT~Time');
% end
% hold on; errorbar(edges(1:end-1)+binWidth/2,nanmean(ContraAll),nanstd(ContraAll)./sqrt(sum(~isnan(ContraAll))),'ro','MarkerSize',10);
% xlim([0,maxTime]); title('Pre Ipsi Trials'); xlabel('Time since last stim (ms)'); ylabel('\Delta RT (ms)');
% tbl = table(RT,Time,'VariableNames',{'RT','Time'});
% fitPreI = fitlm(tbl,'RT~Time');

