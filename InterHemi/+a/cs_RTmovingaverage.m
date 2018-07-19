close all, clear all

load('E:\U\Code\+u\MetaData.mat');  % Once done, instead of loading just have it as a function input?

history = 10;
% avgRTs = zeros(numel(SL),history,2);
% numTrials = zeros(numel(SL),1);
% totAvg = zeros(size(avgRTs,2),size(avgRTs,3));
test1 = [];
test2 = [];
pre1 = [];
pre2 = [];
cond1 = [];
cond2 = [];
post1 = [];
post2 = [];
numTrials = [];
RStim = 0;
LStim = 0;
pre = [];
cond = [];
post = [];


% Reusing AppendReactionTimes code
for i = 1:numel(SL)
    
    cond = char(SL(i).Condition);
    
    if(~isempty(str2num(SL(i).Bad)) || length(cond)<6 || ~strcmp(cond(1:6),'Contra') || isempty(SL(i).trig1) ||str2num(char(SL(i).Date))<20170208)
        continue;
    end
    
    trials = [SL(i).lefttrials(:,2);SL(i).righttrials(:,2)];
    times = [SL(i).rts_l;SL(i).rts_r];
    
    %     index = 1:size(trials,1);
    L = zeros(length(SL(i).lefttrials),1); %just flip these to look at right trials
    R = ones(length(SL(i).righttrials),1);
    tag = [L;R];
    type = zeros(size(tag,1),1);
    data = [trials,times,tag,type];
    data = sortrows(data,1);        % sort by trial
    
    % Some code used to characterize the pattern of trials. looks like
    % there aren't too many (if any) trials where there are three left
    % trials in a row, so just use two. (Richy)
    %     history = 2;                    % number of trials to check
    trialtype = zeros(history,1);
    flag = 1;
    for k = history:size(data,1)
        flag = 1;
        if data(k,3)                    %if it's a contra trial, doesn't matter
            continue;
        end
        for j = 1:history-1
            if(data(k-j,3))
                trialtype(j) = trialtype(j)+1; % contra trial is j trials before
                data(k,4) = j;
                flag = 0;
                break;
            end
        end
        if (flag)
            trialtype(end) = trialtype(end)+1;   % all were left trials
            data(k,4) = history;
        end
    end
    
%     total = [];
%     for j = 1:history
%         total(end+1) = sum(data(:,4)==j);
%     end
%     bar(total)

    window = 100;
    trials = trials(~isnan(times))./1000;
    times = times(~isnan(times));
    
    % cumulate into matrix
%     temp = NaN(history,max(trialtype));
%     Mtemp = nanmean(temp,2);
%     T = Mtemp(~isnan(Mtemp));
%     for k = 1:history
%         temp(k,1:trialtype(k)) = data(data(:,4)==k,2)';
%     end
    

    if(mod(i,3) == 1)
%         figure(1)
%         hold on
%         plot(median,'Color',[.5,.5,.5])
%         scatter(1:10,nanmean(temp,2),[],[.5,.5,.5])
%         pre = horzcat(pre,[trials,times]);
        pre(end+1,1,1:length(trials)) = trials';
        pre(end,2,1:length(times)) = times';
    elseif(mod(i,3) == 2)
%         figure(2)
%         hold on
%         plot(median,'Color',[.5,.5,.5])
%         scatter(1:10,nanmean(temp,2),[],[.5,.5,.5])
%         cond = horzcat(cond,[trials,times]);
        cond(end+1,1,1:length(trials)) = trials';
        cond(end,2,1:length(times)) = times';
    else
%         figure(3)
%         hold on
%         plot(median,'Color',[.5,.5,.5])
%         scatter(1:10,nanmean(temp,2),[],[.5,.5,.5])
%         post = horzcat(post,[trials,times]);
        post(end+1,1,1:length(trials)) = trials';
        post(end,2,1:length(times)) = times';
    end
    
    %     for k = 1:history
    % %         avgRTs(i,k,1) = nanmean(data(data(:,4)==k,2));
    % %         totAvg(k,1) = totAvg(k,1)+avgRTs(i,k,1).*size(SL(i).lefttrials,1);
    % %         avgRTs(i,k,2) = nanstd(data(data(:,4)==k,2));
    % %         totAvg(k,2) = (totAvg(k,1)+avgRTs(i,k,2).^2).*size(SL(i).lefttrials,1);
    %         if k == 1       % contra-ipsi trial
    %             test1 = vertcat(test1,data(data(:,4)==k,2));
    %             if(mod(i,3) == 1)
    %                 pre1 = vertcat(pre1,data(data(:,4)==k,2));
    %             elseif(mod(i,3) == 2)
    %                 cond1 = vertcat(cond1,data(data(:,4)==k,2));
    %             else
    %                 post1 = vertcat(post1,data(data(:,4)==k,2));
    %             end
    %         else            % ipsi-ipsi trial
    %             test2 = vertcat(test2,data(data(:,4)==k,2));
    %             if(mod(i,3) == 1)
    %                 pre2 = vertcat(pre2,data(data(:,4)==k,2));
    %             elseif(mod(i,3) == 2)
    %                 cond2 = vertcat(cond2,data(data(:,4)==k,2));
    %             else
    %                 post2 = vertcat(post2,data(data(:,4)==k,2));
    %             end
    %         end
    %     end
    
    if(stimHemi == 0)
        numTrials(end+1) = sum(~isnan(SL(i).rts_l));
    else
        numTrials(end+1) = sum(~isnan(SL(i).rts_r));
    end
    % End of Richy code
    %     disp(SL(i).expid);
    
%     if(numTrials(end) == 0)
%         disp('ERROR');
%         keyboard;
%     end
    
    if any(SL(i).rts_l<=0) | any(SL(i).rts_r<=0), keyboard; end
    
end

pre(pre == 0) = NaN;
cond(cond == 0) = NaN;
post(post == 0) = NaN;

Mpre = [];
Mcond = [];
Mpost = [];

window = 120; %seconds

for i = 1:size(pre,1)
    figure(1)
    hold on
    m = util.MovingMedian(squeeze(pre(i,:,:)),window);
    Mpre(i,1:length(m)) = m(2,:);
    a = m(1,:);
    a = a(~isnan(a));
    b = m(2,:);
    b = b(~isnan(b));
    plot(a,b,'Color',[.5,.5,.5])
    
    figure(2)
    hold on
    m = util.MovingMedian(squeeze(cond(i,:,:)),window);
    Mcond(i,1:length(m)) = m(2,:); 
    a = m(1,:);
    a = a(~isnan(a));
    b = m(2,:);
    b = b(~isnan(b));
    plot(a,b,'Color',[.5,.5,.5])
    
    figure(3)
    hold on
    m = util.MovingMedian(squeeze(post(i,:,:)),window);
    Mpost(i,1:length(m)) = m(2,:);
    a = m(1,:);
    a = a(~isnan(a));
    b = m(2,:);
    b = b(~isnan(b));
    plot(a,b,'Color',[.5,.5,.5])
end

% Lpre = Mpre == 0;
% Lpre(:,1,1) = 0;
% Lcond = Mcond == 0;
% Lcond(:,1,1) = 0;
% Lpost = Mpost == 0;
% Lpost(:,1,1) = 0;

Mpre(Mpre == 0) = NaN;
Mcond(Mcond == 0) = NaN;
Mpost(Mpost == 0) = NaN;

Medpre = nanmedian(Mpre);
Medcond = nanmedian(Mcond);
Medpost = nanmedian(Mpost);

figure(1)
plot(1:window:window*length(Medpre),Medpre,'r','LineWidth',1.5)
% ylim([140,300])
figure(2)
plot(1:window:window*length(Medcond),Medcond,'r','LineWidth',1.5)
% ylim([140,300])
figure(3)
plot(1:window:window*length(Medpost),Medpost,'r','LineWidth',1.5)
% ylim([140,300])


% pre(pre == 0) = NaN;
% cond(cond == 0) = NaN;
% post(post == 0) = NaN;
% 
% Mpre = nanmedian(pre);
% Mpre = Mpre(Mpre~=0);
% figure(1)
% hold on
% plot(Mpre,'r','Linewidth',1.5)
% 
% Mcond = nanmedian(cond);
% Mcond = Mcond(Mcond~=0);
% figure(2)
% hold on
% plot(Mcond,'r','Linewidth',1.5)
% 
% Mpost = nanmedian(post);
% Mpost = Mpost(Mpost~=0);
% figure(3)
% hold on
% plot(Mpost,'r','Linewidth',1.5)

% 
% figure(4)
% CondDiff = Mpre(1:min(length(Mpre),length(Mcond)))-Mcond(1:min(length(Mpre),length(Mcond)));
% plot(CondDiff)
% 
% figure(5)
% PostDiff = Mpre(1:min(length(Mpre),length(Mpost)))-Mpost(1:min(length(Mpre),length(Mpost)));
% plot(PostDiff)


% Mpre = nanmedian(pre,2);
% Tpre = sum(~isnan(pre),2);
% Spre = nanstd(pre,0,2)./sqrt(Tpre);
% Mcond = nanmedian(cond,2);
% Tcond = sum(~isnan(cond),2);
% Scond = nanstd(cond,0,2)./sqrt(Tcond);
% Mpost = nanmedian(post,2);
% Tpost = sum(~isnan(post),2);
% Spost = nanstd(post,0,2)./sqrt(Tpost);
% 
% figure(1)
% bar(1:10,Mpre,'FaceAlpha',0)
% set(gca,'xticklabel', Tpre) 
% % ylim([190 300])
% 
% figure(2)
% bar(1:10,Mcond,'FaceAlpha',0)
% set(gca,'xticklabel', Tcond) 
% % ylim([190 300])
% 
% figure(3)
% bar(1:10,Mpost,'FaceAlpha',0)
% set(gca,'xticklabel', Tpost) 
% % ylim([190 300])
% 
% figure(4)
% for i = 1:10
%     subplot(2,5,i)
%     bar(1:3,[Mpre(i),Mcond(i),Mpost(i)])
%     set(gca,'xticklabel', [Tpre(i),Tcond(i),Tpost(i)])
%     title(i)
% end
% 
% total_data = [];
% g1 = [];
% g2 = {};
% 
% for i = 1:10
%     for j = 1:length(pre)
%         if(~isnan(pre(i,j)))
%             total_data(end+1) = pre(i,j);
%             g1(end+1) = i;
%             g2(end+1) = {'pre'};
%         end
%     end
% end
% 
% for i = 1:10
%     for j = 1:length(cond)
%         if(~isnan(cond(i,j)))
%             total_data(end+1) = cond(i,j);
%             g1(end+1) = i;
%             g2(end+1) = {'cond'};
%         end
%     end
% end
% 
% for i = 1:10
%     for j = 1:length(post)
%         if(~isnan(post(i,j)))
%             total_data(end+1) = post(i,j);
%             g1(end+1) = i;
%             g2(end+1) = {'post'};
%         end
%     end
% end

