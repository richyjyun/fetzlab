function h = PlotExperiment(SL)
%
% plots a full page of details for a single experiment
%
% arb feb 20 2015
%
% Modified to work with Ubi
% RJY March 16 2017

h = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(h,'visible','off');

if ~exist('notes', 'var'), notes = []; end

subplot(5, 4, 1)
text(0, .95, [SL.Animal 'Experiment: ', SL.Date], 'fontsize', 11, 'fontweight', 'bold')
text(0, .8, ['Experiment Type: ' SL.Condition], 'fontsize', 7)
text(0, .48, ['Stim site(s): ' horzcat(SL.Stim_Loc)], 'fontsize', 7)
text(0, .40, ['N Stim: ' num2str(length(SL.trig1))], 'fontsize', 7);
text(0, .32, ['StimAmp: ' SL.Stim_Amp '.  Delay: ' num2str(SL.Stim_Delay) 'ms'], 'fontsize', 7, 'interpreter', 'none')
if(strcmp(SL.Animal,'Ubi'))
    text(0, .24, ['Stimulator : ' char(SL.Stimulator)], 'fontsize', 7);
end
text(0, .16, ['Normalized Delay: ' num2str(SL.NormDelay)], 'fontsize', 7);
text(0, .08, ['Stim Hemi: ' SL.StimHemi], 'fontsize', 7);
% text(0, 0, ['Notes: ' SL.Notes], 'fontsize', 7);
axis off

window = [0,0.8];
subplot(5, 4, 2) % av accel L
[accs, ts] = a.GetAverageAcc(SL,window);
line(ts, accs(:,1), 'color', 'k', 'linewidth', 1.5), hold on
line(ts, accs(:,2), 'color', 'b', 'linewidth', 1.5, 'linestyle', ':'), hold on
line(ts, accs(:,3), 'color', 'r', 'linewidth', 1.5, 'linestyle', '--'), hold off
title('Average Left Acceleration', 'fontsize', 9)
legend({'pre','cond','post'},'Location','Northeast')
xlabel('trial time (sec)', 'fontsize', 7)
set(gca, 'fontsize', 7)              
xlim(window)

subplot(5, 4, 6) % av accel R
line(ts, accs(:,4), 'color', 'k', 'linewidth', 1.5), hold on
line(ts, accs(:,5), 'color', 'b', 'linewidth', 1.5, 'linestyle', ':'), hold on
line(ts, accs(:,6), 'color', 'r', 'linewidth', 1.5, 'linestyle', '--'), hold off
title('Average Right Acceleration', 'fontsize', 9)
%legend({'left', 'right'}, 'location', 'northwest')
xlabel('trial time (sec)', 'fontsize', 7)
set(gca, 'fontsize', 7)              
xlim(window)

if ~isempty(SL.trig1)
    start = SL.trig1(1);
    finish = SL.trig1(end);
    preinds_l = SL.lefttrials(:,1)<start;
    condinds_l = SL.lefttrials(:,1)>=start & SL.lefttrials(:,1)<=finish;
    postinds_l = SL.lefttrials(:,1)>finish;
    preinds_r = SL.righttrials(:,1)<start;
    condinds_r = SL.righttrials(:,1)>=start & SL.righttrials(:,1)<=finish;
    postinds_r = SL.righttrials(:,1)>finish;
    
end


subplot(5, 4, 3) % cdf trial duration L
if ~isempty(SL.trig1)
    [precdf, prebins] = a.GetCDF(diff(SL.lefttrials(preinds_l,:),1,2), 5);
    [postcdf, postbins] = a.GetCDF(diff(SL.lefttrials(postinds_l,:),1,2), 5);
    line(prebins/1000, precdf, 'color', 'k', 'linewidth', 1.5), hold on
    line(postbins/1000, postcdf, 'color', 'r', 'linewidth', 1.5), hold off
    legend({'pre','post'},'Location','Southeast')
else
    [precdf, prebins] = a.GetCDF(diff(SL.lefttrials,1,2), 5);
    line(prebins/1000, precdf, 'color', 'k', 'linewidth', 1.5)  
end
title('CDF Left Trial Durations', 'fontsize', 9)
xlabel('trial time (sec)', 'fontsize', 7)
set(gca, 'fontsize', 7)              
xlim([0 1])
ylim([0 1])


subplot(5, 4, 7) % cdf trial duration R
if ~isempty(SL.trig1)
    [precdf, prebins] = a.GetCDF(diff(SL.righttrials(preinds_r,:),1,2), 5);
    [postcdf, postbins] = a.GetCDF(diff(SL.righttrials(postinds_r,:),1,2), 5);
    line(prebins/1000, precdf, 'color', 'k', 'linewidth', 1.5), hold on
    line(postbins/1000, postcdf, 'color', 'r', 'linewidth', 1.5), hold off
else
    [precdf, prebins] = a.GetCDF(diff(SL.righttrials,1,2), 5);
    line(prebins/1000, precdf, 'color', 'k', 'linewidth', 1.5)
end
title('CDF Right Trial Durations', 'fontsize', 9)
xlabel('trial time (sec)', 'fontsize', 7)
set(gca, 'fontsize', 7)              
xlim([0 1])
ylim([0 1])

 
subplot(5, 4, 4) % cdf RT L
if ~isempty(SL.trig1)
    [precdf, prebins] = a.GetCDF(SL.rts_l(preinds_l), 5);
    [condcdf, condbins] = a.GetCDF(SL.rts_l(condinds_l), 5);
    [postcdf, postbins] = a.GetCDF(SL.rts_l(postinds_l), 5);
    line(prebins/1000, precdf, 'color', 'k', 'linewidth', 1.5), hold on
    line(condbins/1000, condcdf, 'color', 'b', 'linewidth', 1.5, 'linestyle', ':'), hold on
    line(postbins/1000, postcdf, 'color', 'r', 'linewidth', 1.5, 'linestyle', '--'), hold off
else
     [precdf, prebins] = a.GetCDF(SL.rts_l, 5);
     line(prebins/1000, precdf, 'color', 'k', 'linewidth', 1.5)
end
title('CDF Left RT', 'fontsize', 9)
xlabel('trial time (sec)', 'fontsize', 7)
set(gca, 'fontsize', 7)              
xlim([0 .6])
ylim([0 1])


subplot(5, 4, 8) % cdf RT R
if ~isempty(SL.trig1)
    [precdf, prebins] = a.GetCDF(SL.rts_r(preinds_r), 5);
    [condcdf, condbins] = a.GetCDF(SL.rts_r(condinds_r), 5);
    [postcdf, postbins] = a.GetCDF(SL.rts_r(postinds_r), 5);
    line(prebins/1000, precdf, 'color', 'k', 'linewidth', 1.5), hold on
    line(condbins/1000, condcdf, 'color', 'b', 'linewidth', 1.5, 'linestyle', ':'), hold on
    line(postbins/1000, postcdf, 'color', 'r', 'linewidth', 1.5, 'linestyle', '--'), hold off
else
     [precdf, prebins] = a.GetCDF(SL.rts_r, 5);
     line(prebins/1000, precdf, 'color', 'k', 'linewidth', 1.5)
end
title('CDF Right RT', 'fontsize', 9)
xlabel('trial time (sec)', 'fontsize', 7)
set(gca, 'fontsize', 7)              
xlim([0 .6])
ylim([0 1])


subplot(5, 4, 5) % average acceleration during conditioning
if ~isempty(SL.trig1)
    if length(SL.Condition) >=6 && strcmp(SL.Condition(1:6),'Contra')% if right triggered conditioning session
        [~, ~, lag_stimtarget] = u.FindClosestAfter(SL.righttrials(:,1), SL.trig1+25); % hack for those timed that show up at 0 or just before due to alignment
    elseif length(SL.Condition) >=4 && strcmp(SL.Condition(1:4),'Ipsi')% if left triggered conditioning session
        [~, ~, lag_stimtarget] = u.FindClosestAfter(SL.lefttrials(:,1), SL.trig1+25);
    else
        lag_stimtarget=NaN;
    end
    lag_stimtarget = abs(lag_stimtarget) - 25;
    [accs,ts] = a.GetAverageAccAboutTriggers(SL);
    plot(ts, accs(:,1), 'color', 'k', 'linewidth', 1.5), hold on
    plot(ts, accs(:,2), 'color', 'r', 'linewidth', 1.5)
    legend({'left', 'right'}, 'Location','Northwest')
    %[count, bin] = hist(abs(lag_stimtarget)-25, 0:2:50000);
    % [accs, ts] = ana.GetAverageAcc(SL(1,2));  %Nx4 array like [lefttrials_leftacc, lefttrials_rightacc, righttrials_leftacc, righttrials_rightacc]
    % if any(SL(1,1).stimchannels{1}=='L') % left hemisphere stimulated
    %     stimsite = 'L Hem';
    %     if any(SL(1,2).cond=='L') % left trial triggered
    %         stimtrig = 'Left Trial';
    %         accs = accs(:,2);
    %     elseif any(SL(1,2).cond=='R') % right trial triggered
    %         stimtrig = 'Right Trial';
    %         accs = accs(:,4);
    %     else
    %         keyboard
    %     end
    % elseif any(SL(1,1).stimchannels{1}=='R') | strcmp(SL(1,1).monkey, 'Kato')  % show left accel
    %     stimsite = 'R Hem';
    %     if any(SL(1,2).cond=='L') % left trial triggered
    %         stimtrig = 'Left Trial';
    %         accs = accs(:,1);
    %     elseif any(SL(1,2).cond=='R') % right trial triggered
    %         stimtrig = 'Right Trial';
    %         accs = accs(:,3);
    %     else
    %         keyboard
    %     end
    % end
    
    %bar(bin/1000, count, .7, 'facecolor', 'k'), hold on
    %line(ts, .8*max(count)*accs/max(accs(1:round(SL(1).fs))), 'color', 'r', 'linewidth', 2)
    title(['Stimulus Evoked Twitch.' char(10) 'Effective Lag: ' num2str(nanmedian(lag_stimtarget)) char(177) num2str(nanstd(lag_stimtarget)) 'ms'], 'fontsize', 9)
    xlabel('Time from Stimulation (sec)', 'fontsize', 7)
    set(gca, 'fontsize', 7)
end



subplot(5,4,[9 13]) 
% reaction time bar graph
% rts = vertcat(nanmedian(SL(1,1).rts_l(:)),... % building rt vector
%         nanmedian(SL(1,2).rts_l(:)),...
%         nanmedian(SL(1,3).rts_l(:)),...
%         nanmedian(SL(1,1).rts_r(:)),...
%         nanmedian(SL(1,2).rts_r(:)),...
%         nanmedian(SL(1,3).rts_r(:)));
%     
% bar([1 2 3 5 6 7], rts);
% set(gca, 'fontsize', 7)  
% set(findobj(gca,'Type','text'), 'fontsize',7)
% title('Reaction Time Medians'), ylabel('(ms)')
%% this is the old boxplot that I am temp removing 9/21/15
preleft = repmat('pre_L', sum(preinds_l),1); % build grouping variables
condleft = repmat('con_L', sum(condinds_l),1);
postleft = repmat('pos_L', sum(postinds_l),1);
preright = repmat('pre_R', sum(preinds_r),1);
condright = repmat('con_R', sum(condinds_r),1);
postright = repmat('pos_R', sum(postinds_r),1);

rts = vertcat(SL.rts_l(:),... % building rt vector
        SL.rts_r(:));
    
groups = vertcat(   preleft,...
                    condleft,...
                    postleft,...
                    preright,...
                    condright,...
                    postright);

boxplot(rts, groups, 'notch','on'); %ylim([0 600])
set(gca, 'fontsize', 7)  
set(findobj(gca,'Type','text'), 'fontsize',7)
title('Reaction Times'), ylabel('(ms)')

% Plot RT over time
w = 50;
movL = movmedian(SL.rts_l,w,'omitnan');
movR = movmedian(SL.rts_r,w,'omitnan');

subplot(5,4,10)
plot(movL(preinds_l),'linewidth', 1.5, 'color', 'k'); hold on;
plot(sum(preinds_l)+1:sum(preinds_l)+sum(condinds_l),movL(condinds_l),'linewidth', 1.5, 'color', 'b'); hold on;
plot(sum(preinds_l)+sum(condinds_l)+1:sum(preinds_l)+sum(condinds_l)+sum(postinds_l),movL(postinds_l),'linewidth', 1.5, 'color', 'r'); hold off;
set(gca, 'fontsize', 7)  
title('Left Hand RT')
xlim([1,length(movL)])
ylim([min(movL),max(movL)])

subplot(5,4,14)
plot(movR(preinds_r),'linewidth', 1.5, 'color', 'k'); hold on;
plot(sum(preinds_r)+1:sum(preinds_r)+sum(condinds_r),movR(condinds_r),'linewidth', 1.5, 'color', 'b'); hold on;
plot(sum(preinds_r)+sum(condinds_r)+1:sum(preinds_r)+sum(condinds_r)+sum(postinds_r),movR(postinds_r),'linewidth', 1.5, 'color', 'r'); hold off;
set(gca, 'fontsize', 7)  
title('Right Hand RT')
xlim([1,length(movR)])
ylim([min(movR),max(movR)])

% subplot(5,4,10) % bar plot of betas in both hemispheres
% [preBL, preBR] = ana.CalculateBetaBeforeMovements(SL(1,1), ch_hem=='L'); % betas before left trials, betas before right trials
% [conBL, conBR] = ana.CalculateBetaBeforeMovements(SL(1,2), ch_hem=='L');
% [posBL, posBR] = ana.CalculateBetaBeforeMovements(SL(1,3), ch_hem=='L');
% line(repmat([1;2;3], 1, length(preBL)), [preBL;conBL;posBL], 'color', [.6 .6 .6], 'linewidth', 1, 'marker', 'o'), hold on
% line(repmat([5;6;7], 1, length(preBL)), [preBR;conBR;posBR], 'color', [.6 .6 .6], 'linewidth', 1, 'marker', 'o')
% bar([1 2 3 5 6 7], median([[preBL;conBL;posBL];[preBR;conBR;posBR]],2)', .5, 'edgecolor', 'k', 'facecolor', 'none', 'linewidth', 1.5)
% title('Left Hem \beta'), set(gca, 'xtick', [1 2 3 5 6 7], 'xticklabel', {'PreL', 'ConL', 'PosL', 'PreR', 'ConR', 'PosR'})
% 
% subplot(5,4,14) % bar plot of betas in both hemispheres
% [preBL, preBR] = ana.CalculateBetaBeforeMovements(SL(1,1), ch_hem=='R'); % betas before left trials, betas before right trials
% [conBL, conBR] = ana.CalculateBetaBeforeMovements(SL(1,2), ch_hem=='R');
% [posBL, posBR] = ana.CalculateBetaBeforeMovements(SL(1,3), ch_hem=='R');
% line(repmat([1;2;3], 1, length(preBL)), [preBL;conBL;posBL], 'color', [.6 .6 .6], 'linewidth', 1, 'marker', 'o'), hold on
% line(repmat([5;6;7], 1, length(preBL)), [preBR;conBR;posBR], 'color', [.6 .6 .6], 'linewidth', 1, 'marker', 'o')
% bar([1 2 3 5 6 7], median([[preBL;conBL;posBL];[preBR;conBR;posBR]],2)', .5, 'edgecolor', 'k', 'facecolor', 'none', 'linewidth', 1.5)
% title('Right Hem \beta'), set(gca, 'xtick', [1 2 3 5 6 7], 'xticklabel', {'PreL', 'ConL', 'PosL', 'PreR', 'ConR', 'PosR'})
% 
subplot(5, 4, [11 15]) % histogram of RTs L
bins = 0:8:800;
pre = hist(SL.rts_l(preinds_l), bins);
cond = hist(SL.rts_l(condinds_l), bins) + max(pre);
post = hist(SL.rts_l(postinds_l), bins) + max(cond);

stairs(bins, pre, 'linewidth', 1.5, 'color', 'k'), hold on
stairs(bins, cond, 'linewidth', 1.5, 'color', 'b'), hold on
stairs(bins, post, 'linewidth', 1.5, 'color', 'r'), set(gca, 'xminorgrid', 'on')
set(gca, 'fontsize', 7)  
title('Distribution of RT Left')
% legend('pre','cond','post')

subplot(5, 4, [12 16]) % histogram of RTs R
pre = hist(SL.rts_r(preinds_r), bins);
cond = hist(SL.rts_r(condinds_r), bins) + max(pre);
post = hist(SL.rts_r(postinds_r), bins) + max(cond);

stairs(bins, pre, 'linewidth', 1.5, 'color', 'k'), hold on
stairs(bins, cond, 'linewidth', 1.5, 'color', 'b'), hold on
stairs(bins, post, 'linewidth', 1.5, 'color', 'r'), set(gca, 'xminorgrid', 'on')
set(gca, 'fontsize', 7)  
title('Distribution of RT Right')

preend_l = find(preinds_l,1,'last');
preend_r = find(preinds_r,1,'last');
condend_l = find(condinds_l,1,'last');
condend_r = find(condinds_r,1,'last');
subplot(5, 4, 17:20) % performance
dur_pre = max(SL.lefttrials(preend_l,2), SL.righttrials(preend_r,2));
dur_cond = max(SL.lefttrials(condend_l,2), SL.righttrials(condend_r,2));
tssuccesses = sort([    SL.lefttrials(logical(SL.lefttrialsuccess(preinds_l)),2);...
                        SL.righttrials(logical(SL.righttrialsuccess(preinds_r)),2);...
                        SL.lefttrials(logical(SL.lefttrialsuccess(condinds_l)),2)+dur_pre;...
                        SL.righttrials(logical(SL.righttrialsuccess(condinds_r)),2)+dur_pre;...
                        SL.lefttrials(logical(SL.lefttrialsuccess(postinds_l)),2)+dur_cond;...
                        SL.righttrials(logical(SL.righttrialsuccess(postinds_r)),2)+dur_cond]);

[tfsuccess, bins] = hist(tssuccesses,500000);
window = round(.5*60*1000/diff(bins(1:2))); % 5 minute bin
Ncorrect_t = filter(ones(window,1),1,tfsuccess); % sliding window correct in time
line(bins(bins<dur_pre)/(1000*60), Ncorrect_t(bins<dur_pre), 'color', 'k', 'linewidth', 1.5), hold on
line(bins(bins>dur_pre & bins<dur_cond)/(1000*60), Ncorrect_t(bins>dur_pre & bins<dur_cond), 'color', 'b', 'linewidth', 1.5), hold on
line(bins(bins>dur_cond)/(1000*60), Ncorrect_t(bins>dur_cond), 'color', 'r', 'linewidth', 1.5)
xlabel('Minutes')
ylabel('Correct Trials (/30 sec)')
title('Task Performance During Preconditioning, Conditioning and Postconditioning', 'interpreter', 'none', 'fontsize', 9)
set(gca, 'fontsize', 7)
