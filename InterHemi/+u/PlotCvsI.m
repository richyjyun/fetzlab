% load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbi.mat')
% SL = a.RemoveTwitch(SL); SL = u.AppendThreshold(SL); SL = a.AppendRT(SL); SL = a.AppendNormalizedDelay(SL); UbiSL = SL;
% load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataIgor.mat')
% SL = a.RemoveTwitch(SL); SL = u.AppendThreshold(SL); SL = a.AppendRT(SL); SL = a.AppendNormalizedDelay(SL); IgorSL = SL;
% load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataKato.mat')
% SL = a.RemoveTwitch(SL); SL = u.AppendThreshold(SL); SL = a.AppendRT(SL); SL = a.AppendNormalizedDelay(SL); KatoSL = SL;

figure;
AllRT = cell([2,1]); linfit = struct;
for S = 3
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
        Bad = [20120515,20170518,20170601,20170505];%[20120508,20120515,20120612,20120814];
        if( any(Bad==D))
            continue;
        end
        
        %         Hundred = [20170302,20170303,20170306];
        %         if( ~any(Hundred==D))
        %             continue;
        %         end
        
%         
        if(SL(i).NormDelay < 0)
            continue;
        end
        if(~strcmp(SL(i).Condition(1:6),'Contra'))% || str2num(SL(i).Stim_Delay)==600)
            continue;
        end
        
        disp(SL(i).Date)
        
        %determine contra and ipsi trials
        Contra = []; Ipsi = []; CRT = []; IRT = []; StimTrials = []; Trials = [];
        RT = []; Label = [];
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
        
        
        StimTrials = Contra;
        
        StimInd = zeros(length(StimTrials),1);
        for j = 1:length(SL(i).trig1)
            norm = StimTrials(:,1)-50-SL(i).trig1(j);
            ind = find(norm>0,1)-1;
            if(~isempty(ind))
                StimInd(ind) = 1;
            end
        end  
        
        Trials = [Contra;Ipsi];
        StimInd = [StimInd;zeros(length(Trials)-length(StimInd),1)];
        CRT = CRT - nanmean(CRT(Contra(:,2)<SL(i).trig1(1)));
        IRT = IRT - nanmean(IRT(Ipsi(:,2)<SL(i).trig1(1)));
        RT = [CRT;IRT];
        Label = [ones(length(CRT),1);zeros(length(IRT),1)];  %uint16 code for C and I
            
        %set up order of all trials
        [Trials,Order] = sort(Trials); Order = Order(:,1);
        Trim = [find(Trials(:,2)<SL(i).trig1(1),1,'last'),find(Trials(:,1)>SL(i).trig1(end),1)];
%         Trim = [1,find(Trials(:,2)<SL(i).trig1(1),1,'last')-1];
        RT = RT(Order); RT = RT(Trim(1):Trim(2));
        Label = Label(Order); Label = Label(Trim(1):Trim(2));
        StimInd = StimInd(Order); StimInd = StimInd(Trim(1):Trim(2));
        
        
        StimHist = zeros(length(StimInd),1);
        for j = 2:length(StimHist)
            last = find(StimInd(1:j-1) == 1,1,'last');
            if(~isempty(last))
                StimHist(j) = j-last;
            end
        end
        
        
        
        FindPairs = diff(Label); % difference between consecutive pairs
        FindPairs = (FindPairs==-1);%| FindPairs==-1);  % Can modulate this to compare C to I and I to C.
%         PairInd = find(FindPairs); interval = diff(PairInd);   % interval of the changes, make sure we don't count the same trial multiple times?
%         FindPairs(PairInd(find(interval<2)+1)) = 0; % remove those that aren't unique pairs
        for k = 2:length(FindPairs)
            if(FindPairs(k) && FindPairs(k-1))
                FindPairs(k) = 0;
            end
        end

        FindPairs(end+1) = 0;             % for indexing purposes (needs to be same length as Label).
             
        FindPairs(find(FindPairs)+1) = 1; % all good inds
        
        C = RT(Label == 1 & FindPairs); %CHist = StimHist(Label == 1 & FindPairs);
        I = RT(Label == 0 & FindPairs); %IHist = StimHist(Label == 0 & FindPairs);
        
%         outliers = (abs(C-nanmean(C)) > 3*nanstd(C)) | (abs(I-nanmean(I)) > 3*nanstd(I));
%         C = C(~outliers); I = I(~outliers);
        
%         C = C(1:min(length(C),length(I)));        I = I(1:min(length(C),length(I)));
        
        hold on; scatter(C,I,'k.');
        AllRT{1} = [AllRT{1};C];
        AllRT{2} = [AllRT{2};I];
        tbl = table(C,I,'VariableNames',{'CRT','IRT'});
        linfit(end+1).f = fitlm(tbl,'IRT~CRT');
    end
    
end

xl = xlim; yl = ylim; 
plot(xlim,[0,0],'k--',[0,0],ylim,'k--');
xlabel('Contra \Delta RT (ms)'); ylabel('Ipsi \Delta RT (ms)'); 
tbl = table(AllRT{1},AllRT{2},'VariableNames',{'CRT','IRT'});
fit1 = fitlm(tbl,'IRT~CRT')

Q(1) = sum(AllRT{1}>0 & AllRT{2}>0);
Q(2) = sum(AllRT{1}<0 & AllRT{2}>0);
Q(3) = sum(AllRT{1}<0 & AllRT{2}<0);
Q(4) = sum(AllRT{1}>0 & AllRT{2}<0);
