% load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbi.mat')
% SL = a.RemoveTwitch(SL); SL = u.AppendThreshold(SL); SL = a.AppendRT(SL); SL = a.AppendNormalizedDelay(SL); UbiSL = SL;
% load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataIgor.mat')
% SL = a.RemoveTwitch(SL); SL = u.AppendThreshold(SL); SL = a.AppendRT(SL); SL = a.AppendNormalizedDelay(SL); IgorSL = SL;
% load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataKato.mat')
% SL = a.RemoveTwitch(SL); SL = u.AppendThreshold(SL); SL = a.AppendRT(SL); SL = a.AppendNormalizedDelay(SL); KatoSL = SL;


ContraRT = {};
IpsiRT = {};

NormDelay = [];

for S = 2
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

%         Hundred = [20170302,20170303,20170306];
%         if( ~any(Hundred==D))
%             continue;
%         end
        
        
        if(SL(i).NormDelay < 0)
            continue;
        end
        if(~strcmp(SL(i).Condition(1:6),'Contra'))% || str2num(SL(i).Stim_Delay)==600)
            continue;
        end
        
        NormDelay(end+1) = SL(i).NormDelay;
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
        Trials = [Contra(:,1);Ipsi(:,1)];
        RT = [CRT;IRT];
        Label = char([ones(length(CRT),1)*67;ones(length(IRT),1)*73]);  %uint16 code for C and I
        
        % for conditioning period, only include trials 2 after stim.
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
        Label = Label(Order);
        RT = RT(Order);
        
        %label each trial with how many trials since last stim
        StimHist = zeros(length(StimInd),1);  
        for j = 2:length(StimHist)
            last = find(StimInd(1:j-1) == 1,1,'last');
            if(~isempty(last))
                StimHist(j) = j-last;
            end
        end
        
        %determine cutoff trials
        ContraCut(1) = find(Contra(:,1)>SL(i).trig1(1),1)-1;
        ContraCut(2) = find(Contra(:,1)>SL(i).trig1(end),1);
        IpsiCut(1) = find(Ipsi(:,1)>SL(i).trig1(1),1)-1;
        IpsiCut(2) = find(Ipsi(:,1)>SL(i).trig1(end),1);
        
        %Store RT
        CheckHist = 2;
        ContraRT{1,end+1} = CRT(1:ContraCut(1)-1) - nanmean(CRT(1:ContraCut(1)-1));
%         ContraRT{2,end} = CRT(ContraCut(1):ContraCut(2));
        inds = intersect(strfind(Label','C')',find(StimHist>0 & StimHist<=CheckHist));
        ContraRT{2,end} = RT(inds)- nanmean(CRT(1:ContraCut(1)-1));
        ContraRT{3,end} = CRT(ContraCut(2)+1:end)- nanmean(CRT(1:ContraCut(1)-1));
        IpsiRT{1,end+1} = IRT(1:IpsiCut(1)-1)- nanmean(IRT(1:IpsiCut(1)-1));
%         IpsiRT{2,end} = IRT(IpsiCut(1):IpsiCut(2));
        inds = intersect(strfind(Label','I'),find(StimHist>0 & StimHist<=CheckHist));
        IpsiRT{2,end} = RT(inds)- nanmean(IRT(1:IpsiCut(1)-1));
        IpsiRT{3,end} = IRT(IpsiCut(2)+1:end)- nanmean(IRT(1:IpsiCut(1)-1));
                
    end
end

w = 50;
sizes = cellfun(@length,ContraRT); days = size(ContraRT,2);
Pre = nan(days,max(sizes(1,:))); Cond = nan(days,max(sizes(2,:))); Post = nan(days,max(sizes(3,:)));
for i = 1:days
    Pre(i,1:sizes(1,i)) = ContraRT{1,i};
    Cond(i,1:sizes(2,i)) = ContraRT{2,i};
    Post(i,1:sizes(3,i)) = ContraRT{3,i};
end
% figure; plot(movmean((Pre(:,1:median(sizes(1,:))))',w),'k');
% hold on; plot(movmean((Cond(:,1:median(sizes(2,:))))',w),'b');
% hold on; plot(movmean((Post(:,1:median(sizes(3,:))))',w),'r');

all = [nanmean(Pre(:,1:median(sizes(1,:)))),nanmean(Cond(:,1:median(sizes(2,:)))),nanmean(Post(:,1:median(sizes(3,:))))];
all = movmean(all,w);
l = 1; r = length(nanmean(Pre(:,1:median(sizes(1,:))))); figure; plot(l:r,all(l:r),'k');
l = r+1; r = r+length(nanmean(Cond(:,1:median(sizes(2,:))))); hold on; plot(l:r,all(l:r),'b');
l = r+1; r = r+length(nanmean(Post(:,1:median(sizes(3,:))))); hold on; plot(l:r,all(l:r),'r');


sizes = cellfun(@length,IpsiRT); days = size(IpsiRT,2);
Pre = nan(days,max(sizes(1,:))); Cond = nan(days,max(sizes(2,:))); Post = nan(days,max(sizes(3,:)));
for i = 1:days
    Pre(i,1:sizes(1,i)) = IpsiRT{1,i};
    Cond(i,1:sizes(2,i)) = IpsiRT{2,i};
    Post(i,1:sizes(3,i)) = IpsiRT{3,i};
end
% figure; plot(movmean(nanmean(Pre(:,1:median(sizes(1,:)))),w),'k');
% hold on; plot(movmean(nanmean(Cond(:,1:median(sizes(2,:)))),w),'b');
% hold on; plot(movmean(nanmean(Post(:,1:median(sizes(3,:)))),w),'r');


all = [nanmean(Pre(:,1:median(sizes(1,:)))),nanmean(Cond(:,1:median(sizes(2,:)))),nanmean(Post(:,1:median(sizes(3,:))))];
all = movmean(all,w);
l = 1; r = length(nanmean(Pre(:,1:median(sizes(1,:))))); figure; plot(l:r,all(l:r),'k');
l = r+1; r = r+length(nanmean(Cond(:,1:median(sizes(2,:))))); hold on; plot(l:r,all(l:r),'b');
l = r+1; r = r+length(nanmean(Post(:,1:median(sizes(3,:))))); hold on; plot(l:r,all(l:r),'r');



