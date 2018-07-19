load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbiFinal.mat')
SL = u.AppendThreshold(SL); SL = a.AppendRT(SL); SL = a.AppendNormalizedDelay(SL); 
save('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbiFinal.mat','SL');
load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataIgorFinal.mat')
SL = u.AppendThreshold(SL); SL = a.AppendRT(SL); SL = a.AppendNormalizedDelay(SL); IgorSL = SL;
save('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataIgorFinal.mat','SL');
load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataKatoFinal.mat')
SL = u.AppendThreshold(SL); SL = a.AppendRT(SL); SL = a.AppendNormalizedDelay(SL); KatoSL = SL;
save('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataKatoFinal.mat','SL');


ContraRT = {};
IpsiRT = {};
CondEpoch = [];

NormDelay = [];

interval = 1000; %ms. interval of points in interpolation
w = 120; %in seconds. envelope width
w = w*1000/interval; %converting to s


for S = 1:2
    switch S
        case 1
%             prevInd = 0;
            SL = UbiSL;
        case 2
%             prevInd = length(UbiSL);
            SL = IgorSL;
        case 3
%             prevInd = length(UbiSL)+length(IgorSL);
            SL = KatoSL;
    end
    for i = 1:length(SL)
        if (~isempty(SL(i).Bad) || strcmp(SL(i).Condition, 'Control')...
                || strcmp(SL(i).Condition, 'NaN') || strcmp(SL(i).Condition, 'nostim')...
                || isempty(SL(i).trig1) || strcmp(SL(i).Condition(end),'R'))
            continue;
        end

        if(SL(i).NormDelay < 0)
            continue;
        end
        if(length(SL(i).Condition) < 6 || ~strcmp(SL(i).Condition(1:6),'Contra'))% || str2num(SL(i).Stim_Delay)==600)
            continue;
        end
        
        % For choosing which hemisphere it was stimulated 
        if(~strcmp(SL(i).StimHemi,'R'))
            continue;
        end
        
%         if(~(strcmp(SL(i).Condition,'Control') || strcmp(SL(i).Condition,'NaN') || strcmp(SL(i).Condition,'nostim') ...
%                  || strcmp(SL(i).Condition(end),'R'))) 
%             continue;
%         end
%         if(length(SL(i).accel_raw_l)/SL(i).fs<3500)
%             continue;
%         end

        D = str2num(SL(i).Date);
        Bad = [20120515,20170518];%[20120508,20120515,20120612,20120814];
        if( any(Bad==D))
            continue;
        end
        
        %         Hundred = [20170302,20170303,20170306];
        %         if( ~any(Hundred==D))
        %             continue;
        %         end
      
        
%         NormDelay(end+1) = SL(i).NormDelay;
        disp(SL(i).Date)
        
        %determine contra and ipsi trials
        Contra = []; Ipsi = []; CRT = []; IRT = []; StimTrials = []; Trials = [];
        RT = []; Label = [];
        if(strcmp(SL(i).StimHemi,'L') || strcmp(SL(i).StimHemi,'NaN'))
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
        
        %% Determine trials 1 or 2 away from stim
        %determine which trials have stimulation
        StimTrials = Contra;
        Trials = [Contra(:,1);Ipsi(:,1)];
        RT = [CRT;IRT];
        Label = char([ones(length(CRT),1)*67;ones(length(IRT),1)*73]);  %uint16 code for C and I
        
        % set triggers. if control, set time. Else, just use trig1;
        trig = [];
        if(strcmp(SL(i).Condition,'Control')|| strcmp(SL(i).Condition,'NaN') ||strcmp(SL(i).Condition,'nostim'))
            [Trials,~] = sort(Trials);
     
             trig(1) = Trials(floor(length(Trials)/3),1);
                trig(2) = Trials(floor(2*length(Trials)/3),1);
%             if strcmp(SL(i).Animal, 'Kato')
%                 trig(1) =10*60*SL(i).fs;
%                 trig(2) = 30*60*SL(i).fs;
%             elseif strcmp(SL(i).Animal, 'Igor')
%                 trig(1) = 45*60*1000;
%                 trig(2) = 75*60*1000;
%             else
%                 trig(1) = 25*60*SL(i).fs;
%                 trig(2) = 60*60*SL(i).fs;
%             end
            CKeep = ones(1,length(CRT)); IKeep = ones(1,length(IRT));
        else
            StimInd = zeros(length(StimTrials),1);
            for j = 1:length(SL(i).trig1)
                norm = StimTrials(:,1)-50-SL(i).trig1(j);
                ind = find(norm>0,1)-1;
                if(~isempty(ind))
                    StimInd(ind) = 1;
                end
            end
            StimInd = [StimInd;zeros(length(Trials)-length(StimInd),1)];
            
            maxHist = 2;
            %set up order of all trials
            [Trials,Order] = sort(Trials);
            StimInd = StimInd(Order);
            
            StimHist = ones(length(StimInd),1); Stim = find(StimInd);
            for j = find(StimInd,1):find(StimInd,1,'last')
                last = Stim(find(Stim < j,1,'last'));
                if(~isempty(last))
                    StimHist(j) = j - last;
                end
            end
            Keep = StimHist == 1 | StimHist == 2;
            [~,Return] = sort(Order);
            Keep = Keep(Return);
            CKeep = Keep(1:length(Contra)); IKeep = Keep(length(Contra)+1,end);
            trig = SL(i).trig1;
        end
       
        
        %% Store data
        CRT = CRT - nanmean(CRT(Contra(:,2)<trig(1)));
        IRT = IRT - nanmean(IRT(Ipsi(:,2)<trig(1)));
        if(maxHist > 0)
            CRT(~CKeep) = nan; IRT(~IKeep) = nan;
        end
        % interpolate and smooth
        x = 0:interval*SL(i).fs/1000:Contra(end,1)*SL(i).fs/1000;
        y = interp1(Contra(~isnan(CRT),1)*SL(i).fs/1000,CRT(~isnan(CRT)),x);
        ContraRT{end+1} = movmean(y,w,'omitnan');
        %plot(Contra(:,1)*SL(i).fs/1000,CRT,'o',x,y,':.',x,ContraRT{end},'k');
        x = 0:interval*SL(i).fs/1000:Ipsi(end,1)*SL(i).fs/1000;
        y = interp1(Ipsi(~isnan(IRT),1)*SL(i).fs/1000,IRT(~isnan(IRT)),x);
        IpsiRT{end+1} = movmean(y,w);
        %find left / right cut off of conditioning
        CondEpoch(end+1,:) = floor([trig(1),trig(end)]/interval);
        
    end
end

% figure;
[Contra, Cbounds,CError] = plotRT(ContraRT,CondEpoch,interval,'Contra');
[Ipsi,Ibounds,IError] = plotRT(IpsiRT,CondEpoch,interval,'Ipsi');

function [Line,bounds,Error] = plotRT(RT,Epoch,interval,side)

% if(strcmp(side,'Contra'))
%     subplot(2,2,1);
% else
%     subplot(2,2,2);
% end
Split = {};
t = 0; maxt = 0;
for i = 1:length(RT)
    pre = RT{i}(1:Epoch(i,1)-1); Split{1,end+1} = pre;
    cond = RT{i}(Epoch(i,1):Epoch(i,2));  Split{2,end} = cond;
    post = RT{i}(Epoch(i,2)+1:end);  Split{3,end} = post;
%     disp(num2str(nanmedian(cond)-nanmedian(pre)));
    t = t+nanmedian(cond)-nanmedian(pre);
    %disp([num2str(nanmedian(pre)),' ',num2str(nanmedian(cond)),' ',num2str(nanmedian(post))])
    hold on; l=1; r = length(pre); %plot((l:r)*interval/1000,pre,'k');
    hold on; l=r+1; r = r+length(cond); %plot((l:r)*interval/1000,cond,'b');
    hold on; l=r+1; r = r+length(post); %plot((l:r)*interval/1000,post,'r');
    maxt = max([maxt,r(end)*interval/1000]);
end
% disp(['total ',num2str(t/length(RT))]);
% title([side, ' All Days']); xlabel('Time (s)'); ylabel('\Delta RT');

l = [];
l(1) = min(cellfun(@length,Split(1,:)));
l(2) = min(cellfun(@length,Split(2,:)));
l(3) = min(cellfun(@length,Split(3,:)));
Pre = nan(size(Split,2),l(1));
Cond = nan(size(Split,2),l(2));
Post = nan(size(Split,2),l(3));
for i = 1:size(Split,2)
    Pre(i,1:l(1)) = Split{1,i}(1:l(1));
    Cond(i,1:l(2)) = Split{2,i}(1:l(2));
    Post(i,1:l(3)) = Split{3,i}(1:l(3));
end
PreError = nanstd(Pre);%./sqrt(sum(~isnan(Pre)));
CondError = nanstd(Cond);%./sqrt(sum(~isnan(Cond)));
PostError = nanstd(Post);%./sqrt(sum(~isnan(Post)));
% xlim([0,maxt]); xl = xlim; plot(xl,[0,0],'k--');
Error = [PreError,CondError,PostError];

if(strcmp(side,'Contra'))
    subplot(2,2,3);
else
    subplot(2,2,4);
end
lo = nanmean(Pre)-PreError; hi = nanmean(Pre)+PreError; lo(isnan(lo)) = 0; hi(isnan(hi)) = 0;
l=1; r = length(Pre); p = patch([l:r,r:-1:l,l]*interval/1000,[lo,hi(end:-1:1),lo(1)],'k');
set(p, 'facecolor', [0.8 0.8 0.8], 'edgecolor', 'none');
lo = nanmean(Cond)-CondError; hi = nanmean(Cond)+CondError;
hold on; l=r+1; r = r+size(Cond,2); p = patch([l:r,r:-1:l,l]*interval/1000,[lo,hi(end:-1:1),lo(1)],'b');
set(p, 'facecolor', [0.8 0.8 1], 'edgecolor', 'none');
lo = nanmean(Post)-PostError; hi = nanmean(Post)+PostError;
hold on; l=r+1; r = r+size(Post,2); p = patch([l:r,r:-1:l,l]*interval/1000,[lo,hi(end:-1:1),lo(1)],'r');
set(p, 'facecolor', [1 0.8 0.8], 'edgecolor', 'none');

hold on; l=1; r = size(Pre,2); plot((l:r)*interval/1000,nanmean(Pre),'k','linewidth',1.5);
hold on; l=r+1; r = r+size(Cond,2); plot((l:r)*interval/1000,nanmean(Cond),'b','linewidth',1.5);
hold on; l=r+1; r = r+size(Post,2); plot((l:r)*interval/1000,nanmean(Post),'r','linewidth',1.5);
xlim([0,r(end)*interval/1000]); xl = xlim; plot(xl,[0,0],'k--');

title([side,' Trials']); xlabel('Time (s)'); ylabel('\Delta RT');

Line = [nanmean(Pre),nanmean(Cond),nanmean(Post)];
bounds = [size(Pre,2),size(Pre,2)+size(Cond,2)];
end