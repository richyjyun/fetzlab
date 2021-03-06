function SL = AppendReactionTimes2(SL, dwn, plotit)
%
% For each SL index, adds two fields called rt_l and rt_r that are as long
% as left_trials and right_trials, respectively
%
% Currently uses threshold crossing of 2mV, same as Andrew Richardson
%
% arb 13 april 2015

import OOARB.Utils.*

if ~exist('plotit', 'var'), plotit = 0; end

window = .6; % sec
threshold = 15; %2; % mV

% IgorThreshold = -300; % ms (before end of trial). earlier than this wasnt normal (see detectingmovementslinesfiltfiltrichardson.pdf)
% KatoThreshold = -400;

% removed = [];

% gl = 0;
% gr = 0;
% ll = 0;
% lr = 0;

for i = 1:numel(SL)
%     if(SL(i).Bad)
%         continue;
%     end
    D = SL(i).Date;
    Session = char(D);
    disp(['Reaction Time Session ',Session])
    
    RFilter = u.FilterAccForMovementTimes(SL(i).accel_raw_r, SL(i).fs*dwn, 'richardson');
    LFilter = u.FilterAccForMovementTimes(SL(i).accel_raw_l, SL(i).fs*dwn, 'richardson');
   
    % Plotting to compare blanked out vs original
%     if(~isempty(SL(i).trig1) && any(strmatch('Contra',SL(i).Condition)) && ~strcmp(SL(i).Condition(end),'R'))
%         w = 800;
%         figure(i)
%         id = 0:1:w;
%         trialinds = repmat(floor(SL(i).trig1'.*SL(i).fs./1000), length(id), 1) + repmat(id(:), 1, size(SL(i).trig1,1));
%         Snips = RFilter(trialinds);
%         hold on; plot(mean(Snips')); hold off; title(char(SL(i).Date));
%     end
    
    if isempty(RFilter)
        SL(i).rts_l = nan(size(SL(i).lefttrials,1),1);
        SL(i).rts_r = nan(size(SL(i).righttrials,1),1);
        continue;
    end
    
    LFilter = LFilter./SL(i).Max(1);
    RFilter = RFilter./SL(i).Max(2);
    threshold = 1/3;
    
    inds = -window*SL(i).fs*dwn:1:0; % indeces to look backward
    
    % construct left successes
    leftinds = repmat(SL(i).lefttrials(:,1)', length(inds), 1) + repmat(inds(:), 1, size(SL(i).lefttrials,1));
    % construct right successes
    rightinds = repmat(SL(i).righttrials(:,1)', length(inds), 1) + repmat(inds(:), 1, size(SL(i).righttrials,1));
    
    % round to integers
    leftinds = floor(leftinds);
    rightinds = floor(rightinds);
    
    % remove any windows that extend before recording
    badleft = sum(leftinds<=0,1)>0;
    badright = sum(rightinds<=0,1)>0;
    
    % remove any windows that extend after recording
    badleft = badleft | leftinds(end,:)>length(LFilter);
    badright = badright | rightinds(end,:)>length(LFilter);
    
    leftinds(:,badleft) = 1;
    rightinds(:,badright) = 1;
    
    ts = linspace(-window*SL(i).fs*dwn,0, length(inds));
    
    left_move_inds = u.findfirstineachcolumn(LFilter(leftinds)>threshold); % for each column find first threshold crossing
    right_move_inds = u.findfirstineachcolumn(RFilter(rightinds)>threshold);
    
    %     for j=1:size(leftinds,2)
    %         if(max(LFilter(leftinds(:,j)))<10)
    %             left_move_inds(j) = NaN;
    %         end
    %     end
    %
    %     for j=1:size(rightinds,2)
    %         if(max(RFilter(rightinds(:,j)))<10)
    %             right_move_inds(j) = NaN;
    %         end
    %     end
    
    
    badleft = badleft | isnan(left_move_inds); % logical vector where trial is bad
    badright = badright | isnan(right_move_inds);
    
    ts_left = nan(size(left_move_inds));
    ts_right = nan(size(right_move_inds));
    
    ts_left(~badleft) = ts(left_move_inds(~badleft)); % time relative to each trial end
    ts_right(~badright) = ts(right_move_inds(~badright));
    
    %     if strcmp(SL(i).Animal, 'Igor')
    %         ts_left(ts_left<IgorThreshold) = NaN;
    %         ts_right(ts_right<IgorThreshold) = NaN;
    %     elseif strcmp(SL(i).Animal, 'Kato')
    %         ts_left(ts_left<KatoThreshold) = NaN;
    %         ts_right(ts_right<KatoThreshold) = NaN;
    %     end
    
    SL(i).rts_l = diff(SL(i).lefttrials, 1, 2) + ts_left(:);
    SL(i).rts_r = diff(SL(i).righttrials, 1, 2) + ts_right(:);
    SL(i).rts_l = round(SL(i).rts_l*1000/(SL.fs*dwn),1);
    SL(i).rts_r = round(SL(i).rts_r*1000/(SL.fs*dwn),1);
    
    %     cond = char(SL(i).Condition);
    %     if(~(~isempty(str2num(SL(i).Bad)) || length(cond)<6 || ~strcmp(cond(1:6),'Contra') || isempty(SL(i).trig1) ||str2num(char(SL(i).Date))<20170208))
    %         keyboard;
    %     end
    %
    %     if SL(i).expid == 20120607
    %         keyboard
    %     end
    %     if strcmp(char(SL(i).expid),'20120815')
    %         keyboard
    %     end
    %     keyboard
    
    SL(i).rts_l(SL(i).rts_l>600) = NaN; % remove those over 600ms
    SL(i).rts_r(SL(i).rts_r>600) = NaN;
    SL(i).rts_l(SL(i).rts_l<40) = NaN; % remove those under 40ms
    SL(i).rts_r(SL(i).rts_r<40) = NaN;
    
    %     tl = sum(isnan(SL(i).rts_l));
    %     tr = sum(isnan(SL(i).rts_r));
    %     gl = gl+tl;
    %     gr = gr+tr;
    
    %     ll = ll+sum(isnan(SL(i).rts_l))-tl;
    %     lr = lr+sum(isnan(SL(i).rts_r))-tr;
    
    %     l = sum(isnan(SL(i).rts_l))./length(SL(i).rts_l);
    %     r = sum(isnan(SL(i).rts_r))./length(SL(i).rts_r);
    %     removed(i,1) = l;
    %     removed(i,2) = r;
    %     disp([num2str(l) ' removed for left, ' num2str(r) ' removed for right for session ' num2str(i)]);
    
    if any(SL(i).rts_l<=0) | any(SL(i).rts_r<=0), keyboard; end
    
end

% figure
% histogram(removed(:,1),'Normalization','probability')
% xlabel('Percent Removed')
% ylabel('Percent of Sessions')
% title('Left')
%
% figure
% histogram(removed(:,2),'Normalization','probability')
% xlabel('Percent Removed')
% ylabel('Percent of Sessions')
% title('Right')

if plotit,
    
    packetname = 'posterfigures.ps';
    
    fs = SL(i).fs; % Hz
    
    windowsize = [-.2 10] * fs; % samples, take a wide window to see more
    
    epochs = fs*vertcat(SL(i).lefttrials, SL(i).righttrials)/1000;
    
    acceldata = horzcat(SL(i).accel_raw_l,SL(i).accel_raw_r);
    triggersignals = horzcat(LFilter(:), RFilter(:));
    
    plottingwindows = cat(2, epochs(:,1) + windowsize(1), epochs(:,1) + windowsize(2)); % windows based around target presentation
    
    %normalize to make plotting easier
    acceldata = acceldata-repmat(mean(acceldata,1),size(acceldata,1),1);
    acceldata = acceldata ./ repmat(max(acceldata,[],1), size(acceldata,1),1);
    thresh = threshold ./ max(triggersignals,[],1);
    triggersignals = triggersignals ./ repmat(max(triggersignals,[],1), size(triggersignals,1),1);
    
    epochsinc = floor(size(plottingwindows,1)/60); % sample across the whole recording
    
    figure('paperorientation', 'landscape', 'paperposition', [-1 -.7 12.5 9.6], 'position', [0 0 1100 850])
    
    subplotinds = repmat([1 3 5], 1, 40);
    
    for j = 0:39 % iterate through the plots
        
        samples = plottingwindows(5+j*epochsinc,:); % window
        
        %build the patches
        XS = cat(3, epochs,epochs);
        XS = permute(XS,[3 2 1]);
        XS = reshape(XS, 4, numel(XS)/4,1);
        YS = repmat([-1; 1; 1; -1], 1, size(XS,2));
        
        subaxis(6, 1, subplotinds(j+1), 'padding', .01, 'spacing', .01)
        patch(XS, YS, [.9 .9 .9], 'edgecolor', 'none'), hold on
        line(samples(1):samples(2), acceldata(samples(1):samples(2),1), 'linewidth', 1, 'color', [.4 .4 .4]), hold on
        line(samples(1):samples(2), triggersignals(samples(1):samples(2),1), 'linewidth', 1.5, 'color', [0 0 0]), hold on
        line([samples(1), samples(2)], [thresh(1) thresh(1)], 'color', 'r', 'linestyle', ':')
        %   line(samples(1):samples(2), double(Col9Data(samples(1):samples(2), 7))/4, 'linewidth', 1, 'color', 'b')
        
        XS = fs*SL(i).rts_l/1000;
        XS = [XS, XS]';
        YS = repmat([-1;1], 1, size(XS,2));
        
        line(XS, YS, 'linestyle', ':', 'linewidth', .7, 'color', 'r')
        
        XS = fs*SL(i).triggers/1000;
        XS = [XS, XS]';
        YS = repmat([-1;1], 1, size(XS,2));
        
        line(XS, YS, 'linestyle', ':', 'linewidth', .7, 'color', 'g')
        
        xlim(samples)
        
        
        axis off
        
        % plot the other side now
        %build the patches
        XS = cat(3, epochs,epochs);
        XS = permute(XS,[3 2 1]);
        XS = reshape(XS, 4, numel(XS)/4,1);
        YS = repmat([-1; 1; 1; -1], 1, size(XS,2));
        
        subaxis(6, 1, subplotinds(j+1)+1,'padding', .01, 'spacing', .01)
        patch(XS, YS, [.9 .9 .9], 'edgecolor', 'none'), hold on
        line(samples(1):samples(2), acceldata(samples(1):samples(2),2), 'linewidth', 1, 'color', [.4 .4 .4]), hold on
        line(samples(1):samples(2), triggersignals(samples(1):samples(2),2), 'linewidth', 1.5, 'color', [0 0 0]), hold on
        line([samples(1), samples(2)], [thresh(2) thresh(2)], 'color', 'r', 'linestyle', ':')
        
        XS = fs*SL(i).rts_r/1000;
        XS = [XS, XS]';
        YS = repmat([-1;1], 1, size(XS,2));
        
        line(XS, YS, 'linestyle', ':', 'linewidth', .7, 'color', 'r')
        
        XS = fs*SL(i).triggers/1000;
        XS = [XS, XS]';
        YS = repmat([-1;1], 1, size(XS,2));
        
        line(XS, YS, 'linestyle', ':', 'linewidth', .7, 'color', 'g')
        
        xlim(samples)
        
        
        axis off
        
        if mod(j+1,3)==0
            print(gcf, '-dpsc2', packetname, '-append')
            %print(gcf, '-dpdf', [packetname '.pdf'], '-append')
            clf
        end
        
    end
    
    print(gcf, '-dpsc2', packetname, '-append')
    
end