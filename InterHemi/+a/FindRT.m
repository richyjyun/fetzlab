% Plots reaction time across sessions of one stimulation condition to
% compare pre/cond/post.

% clear; close all;
%
% load('E:\U\Code\+u\MetaData.mat');  % Once done, instead of loading just have it as a function input?
%
% SL = ana.AppendReactionTimes(SL);

clearvars -except SL
R_tot = [];
L_tot = [];
RAccelAvg = [];
LAccelAvg = [];
RAccelFilt = [];
LAccelFilt = [];
rt = [];
rt2 = [];
delay = [];
Sessions = cell(0);

for i = 1:length(SL)
    
    cond = char(SL(i).Condition);
%     
%     % For Contra
%     if(~isempty(str2num(SL(i).Bad)) || length(cond)<6 || ~strcmp(cond(1:6),'Contra') || isempty(SL(i).trig1))
%         continue;
%     end
%     
    % For Ipsi
    if(~isempty(str2num(SL(i).Bad)) || length(cond)<6 || ~strcmp(cond(1:4),'Ipsi') || isempty(SL(i).trig1))
        continue;
    end
    
    D = SL(i).Date;
    S = SL(i).Session_Guger_Train;
    
    Session = [char(D),'_',char(S(2))];
    disp(['Session ',Session])
    Sessions(end+1,1) = cellstr(Session);
    Sessions(end,2) = cellstr(SL(i).Condition);
    % u.trainalign3(Session,2);    % Do this separately with a script, let it run overnight. replacing i16 files (train data)
    
    %     [accdata, trig1, trig2, lefttrials, righttrials, fs] = u.LoadTrain([Session,'.i16']); % load data
    %
    %     data = utils.FilterAccForMovementTimes(double(accdata(:,1)), fs, 'richardson');
    %     ts_moveleft = find(diff(data - repmat(mean(data)+1.5*std(data), size(data,1), 1) > 0)>0)+1;
    %     ts_moveleft = 1000*ts_moveleft/fs; % convert to ms from indices
    %
    %     data = utils.FilterAccForMovementTimes(double(accdata(:,2)), fs, 'richardson');
    %     ts_moveright = find(diff(data - repmat(mean(data)+1.5*std(data), size(data,1), 1) > 0)>0)+1;
    %     ts_moveright = 1000*ts_moveright/fs; % convert to ms from indices
    %
    
    
    %     ts_moveleft = ana.MovementTimes(double(accdata(:,1)), 1000); % detect movements in left
    %     ts_moveright = ana.MovementTimes(double(accdata(:,2)), 1000); % detect movements in right
    %
    %     ts_moveright = utils.FindClosestBefore(righttrials(:,2), ts_moveright); % find right movements that happened just before success of trial
    %     ts_moveleft = utils.FindClosestBefore(lefttrials(:,2), ts_moveleft);
    %
    %     ts_moveleft(isnan(ts_moveleft)) = [];
    %     ts_moveleft(ts_moveleft == 0) = [];
    %     ts_moveright(isnan(ts_moveright)) = [];
    %     ts_moveright(ts_moveright == 0) = [];
    %
    %     rts_right = ts_moveright'-righttrials(:,1);
    %     rts_left = ts_moveleft'-lefttrials(:,1);
    %
    
%     if(strcmp(cond(1:7),'Control'))
%         leftbound = SL(i).lefttrials(floor(length(SL(i).lefttrials)/4),2);
%         rightbound = SL(i).lefttrials(floor(length(SL(i).lefttrials)*3/4),2);
%     else
        leftbound = SL(i).trig1(1);
        rightbound = SL(i).trig1(end);
%     end
    % Need to split up to pre/cond/post here before removing anything
    ind_right = zeros(1,length(SL(i).rts_r));
    ind_left = zeros(1,length(SL(i).rts_l));
    ind_right(SL(i).righttrials(:,1)>leftbound & SL(i).righttrials(:,1)<rightbound) = 1;
    ind_right(SL(i).righttrials(:,1)>rightbound) = 2;
    ind_left(SL(i).lefttrials(:,1)>leftbound & SL(i).lefttrials(:,1)<rightbound) = 1;
    ind_left(SL(i).lefttrials(:,1)>rightbound) = 2;
    
    %     rts_right(rts_right<40 | rts_right>600) = NaN; % remove nonsensical ones
    %     rts_left(rts_left<40 | rts_left>600) = NaN; % remove nonsensical ones
    
    R(1) = nanmedian(SL(i).rts_r(ind_right == 0));
    R(2) = nanmedian(SL(i).rts_r(ind_right == 1));
    R(3) = nanmedian(SL(i).rts_r(ind_right == 2));
    R_tot = vertcat(R_tot,R);
    
    L(1) = nanmedian(SL(i).rts_l(ind_left == 0));
    L(2) = nanmedian(SL(i).rts_l(ind_left == 1));
    L(3) = nanmedian(SL(i).rts_l(ind_left == 2));
    L_tot = vertcat(L_tot,L);
    
    %         figure(1)
    %         hold on
    %         plot(R)
    %
    %         figure(2)
    %         hold on
    %         plot(L)
    
    
    %     SL(i).NormDelay = a.AppendNormalizedDelay(SL(i),'L');
    rt(end+1) = R(2)-R(1);
    rt2(end+1) = L(2)-L(1);
    delay(end+1) = SL(i).NormDelay;
    
    %     RAccel = abs(SL(i).accel_raw_r);
    %     LAccel = abs(SL(i).accel_raw_l);
    %     RSnips = [];
    %     LSnips = [];
    %     window = 600;
    %     for j = 1:length(SL(i).righttrials)
    %         RSnips(j,:) = RAccel(SL(i).righttrials(j,1):SL(i).righttrials(j,1)+window);
    %     end
    %     for j = 1:length(SL(i).lefttrials)
    %         LSnips(j,:) = LAccel(SL(i).lefttrials(j,1):SL(i).lefttrials(j,1)+window);
    %     end
    %
    %     idx = size(RAccelAvg,1) +1;
    %     for j = 1:3
    %         RAccelAvg(idx,j,:) = mean(RSnips(ind_right == (j-1),:));
    %         LAccelAvg(idx,j,:) = mean(LSnips(ind_left == (j-1),:));
    %     end
    %
    %     RFilter = utils.FilterAccForMovementTimes(SL(i).accel_raw_r, SL(i).fs, 'richardson');
    %     LFilter = utils.FilterAccForMovementTimes(SL(i).accel_raw_l, SL(i).fs, 'richardson');
    %     RSnips = [];
    %     LSnips = [];
    %     window = 600;
    %     for j = 1:length(SL(i).righttrials)
    %         RSnips(j,:) = RFilter(SL(i).righttrials(j,1):SL(i).righttrials(j,1)+window);
    %     end
    %     for j = 1:length(SL(i).lefttrials)
    %         LSnips(j,:) = LFilter(SL(i).lefttrials(j,1):SL(i).lefttrials(j,1)+window);
    %     end
    %
    %     idx = size(RAccelFilt,1) +1;
    %     for j = 1:3
    %         RAccelFilt(idx,j,:) = mean(RSnips(ind_right == (j-1),:));
    %         LAccelFilt(idx,j,:) = mean(LSnips(ind_left == (j-1),:));
    %     end
    
    
    %
    %     figure
    %     hold on
    %     cdfplot(SL(i).rts_r(ind_right == 0))
    %     cdfplot(SL(i).rts_r(ind_right == 1))
    %     cdfplot(SL(i).rts_r(ind_right == 2))
    %     legend('Pre','Cond','Post')
    %     title([Session, ' Right RT CDF'])
    %
    %     figure
    %     hold on
    %     cdfplot(SL(i).rts_l(ind_left == 0))
    %     cdfplot(SL(i).rts_l(ind_left == 1))
    %     cdfplot(SL(i).rts_l(ind_left == 2))
    %     legend('Pre','Cond','Post')
    %     title([Session, ' Left RT CDF'])
    %
end

figure(1)
scatter(delay,rt)
title('Right Hand')
xlabel('Normalized Stim Delay')
ylabel('Change in Conditioning RT')

figure(2)
scatter(delay,rt2)
title('Left Hand')
xlabel('Normalized Stim Delay')
ylabel('Change in Conditioning RT')

% RZeroed = [];
% LZeroed = [];
% for i = 1:size(RAccelAvg,1)
%     for j = 1:3
%         RZeroed(i,j,:) = abs(RAccelAvg(i,j,:)-RAccelAvg(i,j,1));
%         LZeroed(i,j,:) = abs(LAccelAvg(i,j,:)-LAccelAvg(i,j,1));
%     end
% end

% day = 3;
%
% figure
% hold on
% plot(squeeze(RZeroed(day,1,:)))
% plot(squeeze(RZeroed(day,2,:)))
% plot(squeeze(RZeroed(day,3,:)))
% legend('Pre','Cond','Post')
%
% [c,l] = xcorr(squeeze(RZeroed(day,1,:)),squeeze(RZeroed(day,2,:)));
% rawlag1 = l(find(c==max(c),1));
% [c,l] = xcorr(squeeze(RZeroed(day,1,:)),squeeze(RZeroed(day,3,:)));
% rawlag2 = l(find(c==max(c),1));
%
% figure
% hold on
% plot(squeeze(RAccelFilt(day,1,:)))
% plot(squeeze(RAccelFilt(day,2,:)))
% plot(squeeze(RAccelFilt(day,3,:)))
% legend('Pre','Cond','Post')
%
% [c,l] = xcorr(squeeze(RAccelFilt(day,1,:)),squeeze(RAccelFilt(day,2,:)));
% filtlag1 = l(find(c==max(c),1));
% [c,l] = xcorr(squeeze(RAccelFilt(day,1,:)),squeeze(RAccelFilt(day,3,:)));
% filtlag2 = l(find(c==max(c),1));

%
% figure(1)
% title('R')
% hold on
% plot(mean(R_tot(1:end,:)),'LineWidth',2)
% legend('1','2','3','4','5','6','7','8')
% figure(2)
% title('L')
% hold on
% plot(mean(L_tot(1:end,:)),'LineWidth',2)
% legend('1','2','3','4','5','6','7','8')

