% function [C,I] = PlotHistory(SL)
% 
% load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbiFinal.mat');
% UbiSL = SL;
% % SL = a.RemoveTwitch(SL); SL = u.AppendThreshold(SL); SL = a.AppendRT(SL); SL = a.AppendNormalizedDelay(SL); UbiSL = SL;
% load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataIgorFinal.mat');
% IgorSL = SL;
% % SL = a.RemoveTwitch(SL); SL = u.AppendThreshold(SL); SL = a.AppendRT(SL); SL = a.AppendNormalizedDelay(SL); IgorSL = SL;
% load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataKatoFinal.mat');
% KatoSL = SL;
% % SL = a.RemoveTwitch(SL); SL = u.AppendThreshold(SL); SL = a.AppendRT(SL); SL = a.AppendNormalizedDelay(SL); KatoSL = SL;


% H = figure;
maxHist = 3;
total = length(UbiSL)+length(IgorSL);%+length(KatoSL);
ContraCount = cell([maxHist,1]);%nan(total,maxHist);
IpsiCount = cell([maxHist,1]);%nan(total,maxHist);
ContraTime = cell([maxHist,1]);%nan(total,maxHist);
IpsiTime = cell([maxHist,1]);%nan(total,maxHist);
ContraHist = cell([maxHist,1]);%nan(total,maxHist);
IpsiHist = cell([maxHist,1]);%nan(total,maxHist);
PreContraCount = cell([maxHist,1]);%nan(total,maxHist);
PreIpsiCount = cell([maxHist,1]);%nan(total,maxHist);
PreContraHist = cell([maxHist,1]);%nan(total,maxHist);
PreIpsiHist = cell([maxHist,1]);%nan(total,maxHist);
PreContra = [];
PreIpsi = [];
ContraExp = zeros(total,1);
IpsiExp = zeros(total,1);

for S = 1
    switch S
        case 1
            prevInd = 0;
            SL = UbiSL;
        case 2
            prevInd = length(UbiSL);
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
%         
%         if(str2num(SL(i).Stim_Delay)<-100)
%             continue;
%         end

%         Good = [20170302,20170303,20170306];
%         if( ~any(Good==D))
%             continue;
%         end
        
        %determine contra and ipsi trials
        Contra = []; Ipsi = []; ContraRT = []; IpsiRT = [];
        if(strcmp(SL(i).StimHemi,'L'))
            Contra = SL(i).righttrials;
            Ipsi = SL(i).lefttrials;
            ContraRT = SL(i).rts_r;
            IpsiRT = SL(i).rts_l;
        elseif(strcmp(SL(i).StimHemi,'R'))
            Contra = SL(i).lefttrials;
            Ipsi = SL(i).righttrials;
            ContraRT = SL(i).rts_l;
            IpsiRT = SL(i).rts_r;
        end
        
%         if(nanmedian(ContraRT(Contra(:,2)<SL(i).trig1(1))) > nanmedian(ContraRT(Contra(:,2)>SL(i).trig1(1) & Contra(:,1)<SL(i).trig1(end))))
%             disp(SL(i).Date);
%             continue;
%         end
        
        %determine which trials have stimulation
        if(length(SL(i).Condition) >=6 && strcmp(SL(i).Condition(1:6),'Contra'))
            StimTrials = Contra;
            Trials = [Contra(:,1);Ipsi(:,1)];
            RT = [ContraRT;IpsiRT];
            Label = char([ones(length(ContraRT),1)*67;ones(length(IpsiRT),1)*73]);  %uint16 code for C and I
%             ContraExp(prevInd+i) = 1;
        else%if(strcmp(SL(i).Condition(1:4),'Ipsi'))
            continue;
            StimTrials = Ipsi;
            Trials = [Ipsi(:,1);Contra(:,1)];
            RT = [IpsiRT;ContraRT];
            Label = char([ones(length(IpsiRT),1)*73;ones(length(ContraRT),1)*67]);  %uint16 code for C and I
            IpsiExp(prevInd+i) = 1;
        end
        
        StimInd = zeros(length(StimTrials),1);
        for j = 1:length(SL(i).trig1)
            norm = StimTrials(:,1)-50-SL(i).trig1(j);
            ind = find(norm>0,1)-1;
            if(~isempty(ind))
                StimInd(ind) = 1;
            end
        end
%         StimInd = StimTrials(:,2)+50 > SL(i).trig1(1) & StimTrials(:,1)-50 < SL(i).trig1(end) ;
        StimInd = [StimInd;zeros(length(Trials)-length(StimInd),1)];
        
        %set up order of all trials
        [Trials,Order] = sort(Trials);
        StimInd = StimInd(Order);
        Trim = [find(StimInd,1),find(StimInd,1,'last')+maxHist];
        StimInd = StimInd(Trim(1):Trim(2));
        Trials = Trials(Trim(1):Trim(2));
        PreRT = RT;
        RT = RT(Order); PreRT = RT(1:Trim(1)-1); RT = RT(Trim(1):Trim(2));
        Label = Label(Order); PreLabel = Label(1:Trim(1)-1); Label = Label(Trim(1):Trim(2));
        ContraPre = nanmedian(PreRT(strfind(PreLabel','C')));
        IpsiPre = nanmedian(PreRT(strfind(PreLabel','I')));
        PreContra(end+1) = nanmedian(ContraPre);%nanmedian(RT(strfind(Label','C')))-ContraPre;
        PreIpsi(end+1) = nanmedian(IpsiPre);%nanmedian(RT(strfind(Label','I')))-IpsiPre;
        
        %Store Pre history and counts. History being since last contra
        %trial.
        StimHist = zeros(length(PreLabel),1); all = strfind(PreLabel','C');
        for j = 2:length(StimHist)
            last = all(find(all<j,1,'last'));
            if(~isempty(last))
                StimHist(j) = j-last;
            end
        end
        for j = 1:maxHist
            inds = intersect(strfind(PreLabel','C'),find(StimHist == j));
            PreContraCount{j} = [PreContraCount{j};length(inds)];
            PreContraHist{j} = [PreContraHist{j};nanmedian(PreRT(inds))];%-ContraPre];
            
            inds = intersect(strfind(PreLabel','I'),find(StimHist == j));
            PreIpsiCount{j} = [PreIpsiCount{j};length(inds)];
            PreIpsiHist{j} = [PreIpsiHist{j};nanmedian(PreRT(inds))];%-IpsiPre];
        end
       
        
        %label each trial with how many trials since last stim
        StimHist = zeros(length(StimInd),1); StimTime = nan(length(StimInd),1); 
        for j = 2:length(StimHist)
%             if(~StimInd(j))
            last = find(StimInd(1:j-1) == 1,1,'last');
            StimHist(j) = j-last;
%             end
%             if(StimInd(j) == 1)
%                 StimHist(j) = 0;
%             else
%                 StimHist(j) = StimHist(j-1)+1;
%             end
            StimTime(j) = Trials(j,1) - Trials(last,1);
        end
        
        %Store history and counts
        for j = 1:maxHist
            if(~strcmp(SL(i).Animal,'Ubi'))
                ContraPre = PreContraHist{j}(end);
                IpsiPre = PreIpsiHist{j}(end);
            end
            
            inds = intersect(strfind(Label','C'),find(StimHist == j));
%             inds = intersect(inds,find(~StimInd)); % only count ones that have stim 
            ContraCount{j} = [ContraCount{j};length(inds)];
            ContraHist{j} = [ContraHist{j};nanmedian(RT(inds))];%-ContraPre];
            ContraTime{j} = [ContraTime{j};nanmedian(StimTime(inds))];
            
            inds = intersect(strfind(Label','I'),find(StimHist == j));
            IpsiCount{j} = [IpsiCount{j};length(inds)];
            IpsiHist{j} = [IpsiHist{j};nanmedian(RT(inds))];%-IpsiPre];
            IpsiTime{j} = [IpsiTime{j};nanmedian(StimTime(inds))];
        end
        
        disp([num2str(length(ContraHist{1})),'. ',SL(i).Date,'-',SL(i).Condition]);
    end
end

for i = 1:maxHist
    ContraHist{i} = ContraHist{i}-PreContra';
    IpsiHist{i} = IpsiHist{i}-PreIpsi';
end

% Plot everything
C = figure;
subplot(2,2,1); bar(1:maxHist,cellfun(@nanmedian,ContraCount));
title('Average Contra Trial Counts'); xlabel('Trials Since Stim'); ylabel('Num Trials');
subplot(2,2,3); Total = cellfun(@isnan,ContraHist,'UniformOutput',false); Total = cellfun(@not,Total,'UniformOutput',false); Total = cellfun(@sum,Total);
errorbar(1:maxHist,cellfun(@nanmedian,ContraHist),cellfun(@nanstd,ContraHist)./sqrt(Total)); x = xlim; xlim([0,x(2)+1])
title('Average Contra Reaction Time'); xlabel('Trials Since Stim'); ylabel('\Delta RT (ms)'); hold on;
% errorbar(0,nanmean(PreContra),std(PreContra)/sqrt(length(PreContra))); hold on;
% scatter([1:maxHist,ContraHist]);
subplot(2,2,2); bar(1:maxHist,cellfun(@nanmedian,IpsiCount));
title('Average Ipsi Trial Counts'); xlabel('Trials Since Stim'); ylabel('Num Trials');
subplot(2,2,4); Total = cellfun(@isnan,IpsiHist,'UniformOutput',false); Total = cellfun(@not,Total,'UniformOutput',false); Total = cellfun(@sum,Total);
errorbar(1:maxHist,cellfun(@nanmedian,IpsiHist),cellfun(@nanstd,IpsiHist)./sqrt(Total));  x = xlim; xlim([0,x(2)+1])
title('Average Ipsi Reaction Time'); xlabel('Trials Since Stim'); ylabel('\Delta RT (ms)'); hold on;
% % errorbar(0,nanmean(PreIpsi),std(PreIpsi)/sqrt(length(PreIpsi)));
% a = axes; t1 = title('Contra Experiments');
% a.Visible = 'off'; t1.Visible = 'on';

% figure; subplot(1,2,1); bar(1:maxHist,cellfun(@nanmedian,ContraTime));
% subplot(1,2,2); bar(1:maxHist,cellfun(@nanmedian,IpsiTime));
% 
% % Plot everything
% C = figure;
% subplot(2,2,1); bar(1:maxHist,cellfun(@nanmean,PreContraCount));
% title('Average Contra Trial Counts'); xlabel('Trials Since Stim'); ylabel('Num Trials');
% subplot(2,2,3); Total = cellfun(@isnan,PreContraHist,'UniformOutput',false); Total = cellfun(@not,Total,'UniformOutput',false); Total = cellfun(@sum,Total);
% errorbar(1:maxHist,cellfun(@nanmean,PreContraHist),cellfun(@nanstd,PreContraHist)./sqrt(Total)); x = xlim; xlim([0,x(2)+1])
% title('Average Contra Reaction Time'); xlabel('Trials Since Stim'); ylabel('\Delta RT (ms)'); hold on;
% subplot(2,2,2); bar(1:maxHist,cellfun(@nanmean,PreIpsiCount));
% title('Average Ipsi Trial Counts'); xlabel('Trials Since Stim'); ylabel('Num Trials');
% subplot(2,2,4); Total = cellfun(@isnan,PreIpsiHist,'UniformOutput',false); Total = cellfun(@not,Total,'UniformOutput',false); Total = cellfun(@sum,Total);
% errorbar(1:maxHist,cellfun(@nanmean,PreIpsiHist),cellfun(@nanstd,PreIpsiHist)./sqrt(Total));  x = xlim; xlim([0,x(2)+1])
% title('Average Ipsi Reaction Time'); xlabel('Trials Since Stim'); ylabel('\Delta RT (ms)'); hold on;
% a = axes; t1 = title('Contra Experiments');
% a.Visible = 'off'; t1.Visible = 'on';

% figure;
% RT = []; grp = [];
% for i = 1:maxHist
%     RT = [RT;ContraHist{i}-nanmedian(PreContraHist{i})];
%     grp = [grp;i*ones(length(ContraHist{i}),1)];
% end
% subplot(1,2,1)
% boxplot(RT,grp); %ylim([-60,150]); 
% hold on; xl = xlim; plot(xl,[0,0],'k--');
% 
% 
% RT = []; grp = [];
% for i = 1:maxHist
%     RT = [RT;IpsiHist{i}-nanmedian(PreIpsiHist{i})];
%     grp = [grp;i*ones(length(IpsiHist{i}),1)];
% end
% subplot(1,2,2); 
% boxplot(RT,grp); %ylim([-60,150]);
% hold on; xl = xlim; plot(xl,[0,0],'k--');

% 
% figure;
% scatter([1,2,3],cellfun(@nanmedian,IpsiHist)-cellfun(@nanmedian,PreIpsiHist))
% 

% figure; hold on; all = [];
% for i = 1:maxHist
%     val = ContraHist{i};
%     plot(sort(val))
%     all(i,:) = val;
% end
% legend({'1','2','3','4','5','6'})
% figure; plot(all); 

% I = figure;
% subplot(2,2,1); bar(1:maxHist,nanmean(ContraCount(find(IpsiExp),:)));
% title('Average Contra Trial Counts'); xlabel('Trials Since Stim'); ylabel('Num Trials');
% subplot(2,2,3); y = ContraHist(find(IpsiExp),:);
% errorbar(1:maxHist,nanmean(y),nanstd(y)./sqrt(sum(~isnan(y))));  x = xlim; xlim([x(1)-1,x(2)+1])
% title('Average Contra Reaction Time'); xlabel('Trials Since Stim'); ylabel('RT (ms)');
% subplot(2,2,2); bar(1:maxHist,nanmean(IpsiCount(find(IpsiExp),:)));
% title('Average Ipsi Trial Counts'); xlabel('Trials Since Stim'); ylabel('Num Trials');
% subplot(2,2,4); y = IpsiHist(find(IpsiExp),:);
% errorbar(1:maxHist,nanmean(y),nanstd(y)./sqrt(sum(~isnan(y))));  x = xlim; xlim([x(1)-1,x(2)+1])
% title('Average Ipsi Reaction Time'); xlabel('Trials Since Stim'); ylabel('RT (ms)');
% a = axes; t1 = title('Ipsi Experiments');
% a.Visible = 'off'; t1.Visible = 'on';
    
% end