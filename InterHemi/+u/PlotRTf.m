% Plots reaction time across sessions of one stimulation condition to
% compare pre/cond/post.

% clear; close all;
%
% load('E:\U\Code\+u\MetaData.mat');  % Once done, instead of loading just have it as a function input?
%
% SL = ana.AppendReactionTimes(SL);
function [c,p,g,l,c2,CC_tot,CI_tot] = PlotRTf(SL, Split, c2,g,CC_tot,CI_tot)
%
% SL, split? c2 = figure handle for control fig
%
if(nargin<2)
    Split = 0;
end

stimprob = 0;

switch SL(1).Animal
    case 'Igor'
        constring = 'nostim';
        cs = [.8,.5,.5];
        achar = '*';
        anshift = 0;
    case 'Kato'
        constring = 'NaN';
        cs = [.5,.5,.8];
        achar = '.';
        anshift = -0.1;
    case 'Ubi'
        constring = 'Control';
        cs = [.5,.8,.5];
        achar = '+';
        anshift = .1;
end

% make xs for bar alts
xs = repmat([1 2 3], 2, 1) + repmat([-.2; .2], 1,3) + anshift;

c = figure;
p = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 8]);
if ~exist('g')
    g = figure;
end
l = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 8]);

if ~exist('c2', 'var')
    c2 = figure;
end

if ~exist('CC_tot','var')
    CC_tot = [];
end
CC_tot2 = [];
if ~exist('CI_tot','var')
    CI_tot = [];
end
CI_tot2 = [];
IC_tot = [];
IC_tot2 = [];
II_tot = [];
II_tot2 = [];
CONR_tot = [];
CONL_tot = [];
% rt = [];
% rt2 = [];
% delay = [];
Sessions = cell(0);

for i = 1:length(SL)
    
    cond = char(SL(i).Condition);
    
%     if(~isempty(SL(i).Bad) || isempty(SL(i).trig1)|| strcmp(SL(i).Condition,'nostim') || strcmp(SL(i).Condition,'tonic')...
%             || strcmp(SL(i).Condition,'MOVE_L') || strcmp(SL(i).Condition,'MOVE_R') || strcmp(SL(i).StimHemi,'NaN') ...
%             || strcmp(SL(i).Condition(end),'R'))
% %             || (isfield(SL,'Stimulator') && strcmp(SL(i).Stimulator,'badmcs')))
%         continue;
%     end
     
%     if((str2num(SL(i).Date) >= 20170224 && str2num(SL(i).Date) <= 20170306))
%         continue;
%     end
    if (~isempty(SL(i).Bad) && ~strcmp(SL(i).Condition, 'Control')), continue; end
    
    if(length(cond)>=6 && strcmp(cond(1:6),'Contra') && cond(end)~='R' && ~isempty(SL(i).trig1))
        cond = 0;
    elseif(length(cond)>=4 && strcmp(cond(1:4),'Ipsi') && cond(end)~='R' && ~isempty(SL(i).trig1))
        cond = 1;
    else
        cond = -1;
        if( strcmp(cond, 'nostim') || ( strcmp(cond, 'NaN') && strcmp(SL(1).Animal,'Kato')) )
            achar = '*';
        elseif ( strcmp(cond, 'tonic'))
            achar = '.';
        elseif( strcmp(cond(end),'R'))
            achar = '+';
        else
            continue;
        end
    end
        
    
%         stimSite = 'LMC 4/4.5';
%     if(~strcmp(SL(i).Stim_Loc,stimSite))
%         continue;
%     end
    
    
    D = SL(i).Date;
    
    Session = char(D);
    disp(['Session ',Session])
    Sessions(end+1,1) = cellstr(Session);
    Sessions(end,2) = cellstr(SL(i).Condition);
    if cond==-1 % control, so make up conditionining period
        if strcmp(SL(i).Animal, 'Kato')
            leftbound =10*60*SL(i).fs;
            rightbound = 30*60*SL(i).fs;
        elseif strcmp(SL(i).Animal, 'Igor')
            leftbound = 45*60*1000;
            rightbound = 75*60*1000;
        else
            leftbound = 25*60*SL(i).fs;
            rightbound = 60*60*SL(i).fs;
        end
    elseif ~isempty(SL(i).trig1) % conditioning bounds, I guess
        leftbound = SL(i).trig1(1);
        rightbound = SL(i).trig1(end);
    else
        continue
    end
    
    if cond==-1
        if rand(1)<.5 % randomize whether a control was left or right stim
            IpsiTrials = SL(i).lefttrials; 
            ContTrials = SL(i).righttrials;
            IpsiRT = SL(i).rts_l;
            ContRT = SL(i).rts_r;
        else
            IpsiTrials = SL(i).righttrials; 
            ContTrials = SL(i).lefttrials;
            IpsiRT = SL(i).rts_r;
            ContRT = SL(i).rts_l;
        end
    elseif(strcmp(SL(i).StimHemi,'R'))
        IpsiTrials = SL(i).righttrials;
        ContTrials = SL(i).lefttrials;
        IpsiRT = SL(i).rts_r;
        ContRT = SL(i).rts_l;
    elseif(strcmp(SL(i).StimHemi,'L'))
        IpsiTrials = SL(i).lefttrials;
        ContTrials = SL(i).righttrials;
        IpsiRT = SL(i).rts_l;
        ContRT = SL(i).rts_r;            
    end
         
    ind_Cont = zeros(1,length(ContRT)); % these are inds where contralateral trials are in conditioning
    ind_Ipsi = zeros(1,length(IpsiRT)); % these are inds where ipis...
    ind_Cont(ContTrials(:,1)>leftbound & ContTrials(:,1)<rightbound) = 1;
    ind_Cont(ContTrials(:,1)>rightbound) = 2;
    ind_Ipsi(IpsiTrials(:,1)>leftbound & IpsiTrials(:,1)<rightbound) = 1;
    ind_Ipsi(IpsiTrials(:,1)>rightbound) = 2;
    
    if(strcmp(SL(i).Animal,'Ubi') && cond == 0 && Split) % Contra stim
        stimTrials = zeros(1,length(ContTrials));
        for j = 1:length(ContTrials)
            if(any(SL(i).trig1>ContTrials(j,1)-50 & SL(i).trig1<ContTrials(j,2)+50))
                stimTrials(j) = 1;
            end
        end
        stimprob = sum(stimTrials==1)/sum(ind_Cont==1);
        %disp(['percent stim = ' num2str(stimprob)]);
        R(2) = nanmedian(ContRT(ind_Cont == 1 & ~stimTrials));
        RStim = nanmedian(ContRT(ind_Cont == 1 & stimTrials));
    else
        R(2) = nanmedian(ContRT(ind_Cont == 1));
    end
    R(1) = nanmedian(ContRT(ind_Cont == 0));
    R(3) = nanmedian(ContRT(ind_Cont == 2));
    R = R-50;

    if(strcmp(SL(i).Animal,'Ubi') && cond == 1 && Split) % Ipsi stim
        stimTrials = zeros(1,length(IpsiTrials));
        for j = 1:length(IpsiTrials)
            if(any(SL(i).trig1>IpsiTrials(j,1)-50 & SL(i).trig1<IpsiTrials(j,2)+50))
                stimTrials(j) = 1;
            end
        end
        stimprob = sum(stimTrials==1)/sum(ind_Ipsi==1);
        %disp(['percent stim = ' num2str(stimprob)]);
        L(2) = nanmedian(IpsiRT(ind_Ipsi == 1 & ~stimTrials));
        LStim = nanmedian(IpsiRT(ind_Ipsi == 1 & stimTrials));
    else
        L(2) = nanmedian(IpsiRT(ind_Ipsi == 1));
    end
    L(1) = nanmedian(IpsiRT(ind_Ipsi == 0));
    L(3) = nanmedian(IpsiRT(ind_Ipsi == 2));
    L = L-50;

%     rt(end+1) = R(2)-R(1);
%     rt2(end+1) = L(2)-L(1);
%     delay(end+1) = SL(i).NormDelay;
    
    % Contra exp plots
    if(cond == 0 & ~(Split & stimprob>.8))
        figure(c) % pre - cond
        subplot(1,2,1); hold on;
        scatter(SL(i).NormDelay,R(2)-R(1),'k');
        if ~Split
            err = sqrt( nanstd(ContRT(ind_Cont == 0))^2 / sum(~isnan(ContRT(ind_Cont == 0))) + nanstd(ContRT(ind_Cont == 1))^2 / sum(~isnan(ContRT(ind_Cont == 1))) ); % std of diff the diffmeans (not std error)
        else
            err = sqrt( nanstd(ContRT(ind_Cont == 0))^2 / sum(~isnan(ContRT(ind_Cont == 0))) + nanstd(ContRT(ind_Cont == 1 & ~stimTrials))^2 / sum(~isnan(ContRT(ind_Cont == 1 & ~stimTrials))) );
            err2 = sqrt( nanstd(ContRT(ind_Cont == 0))^2 / sum(~isnan(ContRT(ind_Cont == 0))) + nanstd(ContRT(ind_Cont == 1 & stimTrials))^2 / sum(~isnan(ContRT(ind_Cont == 1 & stimTrials))) );
        end
        errorbar(SL(i).NormDelay,R(2)-R(1),err,'k');
        errorbar(SL(i).NormDelay,R(2)-R(1),SL(i).NormSE,'horizontal','k'); hold off;
%         text(SL(i).NormDelay+10,R(2)-R(1),char(D),'FontSize',7);
        if(strcmp(SL(i).Animal,'Ubi') && Split)
            hold on;
            scatter(SL(i).NormDelay,RStim-R(1),'r')
            errorbar(SL(i).NormDelay,RStim-R(1),err2,'r');
            text(SL(i).NormDelay+10,RStim-R(1),char(D),'FontSize',7,'Color','r');
        end
        
        subplot(1,2,2); hold on;
        scatter(SL(i).NormDelay,L(2)-L(1),'k');
        err = sqrt( nanstd(IpsiRT(ind_Ipsi == 0))^2 / sum(~isnan(IpsiRT(ind_Ipsi == 0))) + nanstd(IpsiRT(ind_Ipsi == 1))^2 / sum(~isnan(IpsiRT(ind_Ipsi == 1))) );
        errorbar(SL(i).NormDelay,L(2)-L(1),err,'k');
        errorbar(SL(i).NormDelay,L(2)-L(1),SL(i).NormSE,'horizontal','k'); hold off;
%         text(SL(i).NormDelay+10,L(2)-L(1),char(D),'FontSize',7);
        
        figure(p) % pre-post
        subplot(2,2,1); hold on;
        scatter(SL(i).NormDelay,R(3)-R(1),'k');
        err = sqrt( nanstd(ContRT(ind_Cont == 0))^2 / sum(~isnan(ContRT(ind_Cont == 0))) + nanstd(ContRT(ind_Cont == 2))^2 / sum(~isnan(ContRT(ind_Cont == 2))) );
        errorbar(SL(i).NormDelay,R(3)-R(1),err,'k');
        errorbar(SL(i).NormDelay,R(3)-R(1),SL(i).NormSE,'horizontal','k'); hold off;
%         text(SL(i).NormDelay+10,R(3)-R(1),char(D),'FontSize',7);
        
        
        subplot(2,2,2); hold on;
        scatter(SL(i).NormDelay,L(3)-L(1),'k');
        err = sqrt( nanstd(IpsiRT(ind_Ipsi == 0))^2 / sum(~isnan(IpsiRT(ind_Ipsi == 0))) + nanstd(IpsiRT(ind_Ipsi == 2))^2 / sum(~isnan(IpsiRT(ind_Ipsi == 2))) );
        errorbar(SL(i).NormDelay,L(3)-L(1),err,'k');
        errorbar(SL(i).NormDelay,L(3)-L(1),SL(i).NormSE,'horizontal','k'); hold off;
%         text(SL(i).NormDelay+10,L(3)-L(1),char(D),'FontSize',7);
        
        if(SL(i).NormDelay>0)
            CC_tot = vertcat(CC_tot,R);
            CI_tot = vertcat(CI_tot,L);
            figure(g)
            subplot(2,2,1); hold on;
            x = [1,2,3];
            plot(x,R,'Color',cs); scatter(x,R,'MarkerEdgeColor',cs);
            subplot(2,2,2); hold on;
            plot(x,L,'Color',cs); scatter(x,L,'MarkerEdgeColor',cs);
        else
            CC_tot2 = vertcat(CC_tot2,R);
            CI_tot2 = vertcat(CI_tot2,L);
            figure(l)
            subplot(2,2,1); hold on;
            x = [1,2,3];
            plot(x,R,'Color',cs); scatter(x,R,'MarkerEdgeColor',cs);
            subplot(2,2,2); hold on;
            plot(x,L,'Color',cs); scatter(x,L,'MarkerEdgeColor',cs);
        end
          
        % Ipsi exp plots
    elseif(cond == 1 & ~(Split & stimprob>.8))
       
        figure(p)
        subplot(2,2,3); hold on;
        scatter(SL(i).NormDelay,R(3)-R(1),'k');
        err = sqrt( nanstd(ContRT(ind_Cont == 0))^2 / sum(~isnan(ContRT(ind_Cont == 0))) + nanstd(ContRT(ind_Cont == 2))^2 / sum(~isnan(ContRT(ind_Cont == 2))) );
        errorbar(SL(i).NormDelay,R(3)-R(1),err,'k');
        errorbar(SL(i).NormDelay,R(3)-R(1),SL(i).NormSE,'horizontal','k'); hold off;
%         text(SL(i).NormDelay+10,R(3)-R(1),char(D),'FontSize',7);
        
        subplot(2,2,4); hold on;
        scatter(SL(i).NormDelay,L(3)-L(1),'k');
        err = sqrt( nanstd(IpsiRT(ind_Ipsi == 0))^2 / sum(~isnan(IpsiRT(ind_Ipsi == 0))) + nanstd(IpsiRT(ind_Ipsi == 2))^2 / sum(~isnan(IpsiRT(ind_Ipsi == 2))) );
        errorbar(SL(i).NormDelay,L(3)-L(1),err,'k');
        errorbar(SL(i).NormDelay,L(3)-L(1),SL(i).NormSE,'horizontal','k'); hold off;
%         text(SL(i).NormDelay+10,L(3)-L(1),char(D),'FontSize',7);
                
        if(SL(i).NormDelay>0)
            IC_tot = vertcat(IC_tot,R);
            II_tot = vertcat(II_tot,L);
            figure(g)
            subplot(2,2,3); hold on;
            x = [1,2,3];
            plot(x,R,'Color',[.5,.5,.5]); scatter(x,R,'MarkerEdgeColor',[.5,.5,.5]);
            subplot(2,2,4); hold on;
            plot(x,L,'Color',[.5,.5,.5]); scatter(x,L,'MarkerEdgeColor',[.5,.5,.5]);
       else
            IC_tot2 = vertcat(IC_tot2,R);
            II_tot2 = vertcat(II_tot2,L);
            figure(l)
            subplot(2,2,3); hold on;
            x = [1,2,3];
            plot(x,R,'Color',[.5,.5,.5]); scatter(x,R,'MarkerEdgeColor',[.5,.5,.5]);
            subplot(2,2,4); hold on;
            plot(x,L,'Color',[.5,.5,.5]); scatter(x,L,'MarkerEdgeColor',[.5,.5,.5]);
        end
        
    elseif cond==-1
        if strcmp(SL(i).Date, '20170405'), continue; end
        CONR_tot = vertcat(CONR_tot,R);
        CONL_tot = vertcat(CONL_tot,L);
        
        figure(c2)
        subplot(2,2,2); hold on;
        x = [1,2,3]+anshift;
        plot(x,R,'Color',cs, 'linewidth', 1.5); scatter(x,R,300,'MarkerEdgeColor',cs, 'marker', achar, 'linewidth', 2);
        subplot(2,2,1); hold on;
        plot(x,L,'Color',cs, 'linewidth', 1.5); scatter(x,L,300,'MarkerEdgeColor',cs, 'marker', achar, 'linewidth', 2);
       
    else
        continue;
    end
    
end


% a = axes; t1 = title('Cond vs Pre','Color','r');
a.Visible = 'off'; t1.Visible = 'on';

figure(p)
subplot(2,2,1)
title('Contra Cond - Contra Trial')
xlabel('Normalized Stim Delay')
ylabel('Change in Conditioning RT')

subplot(2,2,2)
title('Contra Cond - Ipsi Trial')
xlabel('Normalized Stim Delay')
ylabel('Change in Conditioning RT')

subplot(2,2,3)
title('Ipsi Cond - Contra Trial')
xlabel('Normalized Stim Delay')
ylabel('Change in Conditioning RT')

subplot(2,2,4)
title('Ipsi Cond - Ipsi Trial')
xlabel('Normalized Stim Delay')
ylabel('Change in Conditioning RT')

a = axes; t1 = title('Post vs Pre','Color','r');
a.Visible = 'off'; t1.Visible = 'on';

%% bar plots
ylims = vertcat([min(min(CC_tot))-10,max(max(CC_tot))+10],...
                [min(min(CI_tot))-10,max(max(CI_tot))+10],...
                [min(min(IC_tot))-10,max(max(IC_tot))+10],...
                [min(min(II_tot))-10,max(max(II_tot))+10]);
ylims = [min(ylims(:,1)), max(ylims(:,2))];
ylims = [100 ylims(2)]; % set min to 100ms

figure(g)
if(strcmp(SL(1).Animal,'Kato'))
    subplot(2,2,1)
    title('Contra Cond - Contra Trial');
    ylabel('RT (ms)')
    if(~isempty(CC_tot))
        ylim(ylims); hold on;
        b = bar([1,2,3],mean(CC_tot,1),'w');
        uistack(b,'bottom');
    end
    set(gca,'XTick',1:3);
    set(gca,'XTickLabel',{'PRE','STIM','POST'});
    
    subplot(2,2,2)
    title('Contra Cond - Ipsi Trial')
    ylabel('RT (ms)')
    if(~isempty(CI_tot))
        ylim(ylims)
        b = bar([1,2,3],mean(CI_tot,1),'w');
        uistack(b,'bottom');
    end
    set(gca,'XTick',1:3);
    set(gca,'XTickLabel',{'PRE','STIM','POST'});
end

subplot(2,2,3)
title('Ipsi Cond - Contra Trial')
ylabel('RT (ms)')
if(~isempty(IC_tot))
    ylim(ylims)
    b = bar([1,2,3],median(IC_tot,1),'w');
    uistack(b,'bottom');
end
set(gca,'XTick',1:3);
set(gca,'XTickLabel',{'PRE','STIM','POST'});

subplot(2,2,4)
title('Ipsi Cond - Ipsi Trial')
ylabel('RT (ms)')
if(~isempty(II_tot))
    ylim(ylims)
    b = bar([1,2,3],median(II_tot,1),'w');
    uistack(b,'bottom');
end
set(gca,'XTick',1:3);
set(gca,'XTickLabel',{'PRE','STIM','POST'});


figure(l)
subplot(2,2,1)
title('Contra Cond - Contra Trial');
ylabel('RT (ms)')
if(~isempty(CC_tot2))
    ylim(ylims)
    b = bar([1,2,3],median(CC_tot2,1),'w');
    uistack(b,'bottom');
end
set(gca,'XTick',1:3);
set(gca,'XTickLabel',{'PRE','STIM','POST'});

subplot(2,2,2)
title('Contra Cond - Ipsi Trial')
ylabel('RT (ms)')
if(~isempty(CI_tot2))
    ylim(ylims)
    b = bar([1,2,3],median(CI_tot2,1),'w');
    uistack(b,'bottom');
end
set(gca,'XTick',1:3);
set(gca,'XTickLabel',{'PRE','STIM','POST'});

subplot(2,2,3)
title('Ipsi Cond - Contra Trial')
ylabel('RT (ms)')
if(~isempty(IC_tot2))
    ylim(ylims)
    b = bar([1,2,3],median(IC_tot2,1),'w');
    uistack(b,'bottom');
end
set(gca,'XTick',1:3);
set(gca,'XTickLabel',{'PRE','STIM','POST'});

subplot(2,2,4)
title('Ipsi Cond - Ipsi Trial')
ylabel('RT (ms)')
if(~isempty(II_tot2))
    ylim(ylims)
    b = bar([1,2,3],median(II_tot2,1),'w');
    uistack(b,'bottom');
end
set(gca,'XTick',1:3);
set(gca,'XTickLabel',{'PRE','STIM','POST'});

figure(c2)
jitterAmount = 0.5;
jitX = 2*(rand(1,3)-0.5)*jitterAmount;   % +/-jitterAmount max

subplot(2,2,1)
title('CONTROL CONTRA');
ylabel('RT (ms)')
if(~isempty(CONL_tot))
    ylim(ylims)

    plot(xs, repmat(mean(CONL_tot,1),2,1), 'color', .75*cs, 'linewidth', 1.5) 

end
set(gca,'XTick',1:3);
set(gca,'XTickLabel',{'PRE','STIM','POST'});

subplot(2,2,2)
title('CONTROL IPSI')
ylabel('RT (ms)')
if(~isempty(CONR_tot))
    ylim(ylims)

    plot(xs, repmat(mean(CONR_tot,1),2,1), 'color', .75*cs, 'linewidth', 1.5)

end
set(gca,'XTick',1:3);
set(gca,'XTickLabel',{'PRE','STIM','POST'});

% 
% subplot(2,2,1)
% scatter(delay,rt)
% title('Contra Cond - Contra Trial')
% xlabel('Normalized Stim Delay')
% ylabel('Change in Conditioning RT')
% 
% subplot(2,2,2)
% scatter(delay,rt2)
% title('Contra Cond - Ipsi Trial')
% xlabel('Normalized Stim Delay')
% ylabel('Change in Conditioning RT')
% 
% 
% 
% 
% R_tot = [];
% L_tot = [];
% rt = [];
% rt2 = [];
% delay = [];
% Sessions = cell(0);
% for i = 1:length(SL)
%     
%     cond = char(SL(i).Condition);
%     
%     % For Ipsi
%     if(~isempty(SL(i).Bad) || length(cond)<4 || ~strcmp(cond(1:4),'Ipsi') || isempty(SL(i).trig1))
%         continue;
%     end 
%     
%     D = SL(i).Date;
%     
%     Session = char(D);
%     disp(['Session ',Session])
%     Sessions(end+1,1) = cellstr(Session);
%     Sessions(end,2) = cellstr(SL(i).Condition);
%     
%     leftbound = SL(i).trig1(1);
%     rightbound = SL(i).trig1(end);
%     
%     if(strcmp(SL(i).StimHemi,'R'))
%         IpsiTrials = SL(i).righttrials;
%         ContTrials = SL(i).lefttrials;
%         IpsiRT = SL(i).rts_r;
%         ContRT = SL(i).rts_l;
%     elseif(strcmp(SL(i).StimHemi,'L'))
%         IpsiTrials = SL(i).lefttrials;
%         ContTrials = SL(i).righttrials;
%         IpsiRT = SL(i).rts_l;
%         ContRT = SL(i).rts_r;
%     else
%         continue;
%     end
%          
%     ind_Cont = zeros(1,length(ContRT));
%     ind_Ipsi = zeros(1,length(IpsiRT));
%     ind_Cont(ContTrials(:,1)>leftbound & ContTrials(:,1)<rightbound) = 1;
%     ind_Cont(ContTrials(:,1)>rightbound) = 2;
%     ind_Ipsi(IpsiTrials(:,1)>leftbound & IpsiTrials(:,1)<rightbound) = 1;
%     ind_Ipsi(IpsiTrials(:,1)>rightbound) = 2;
%     
%     R(1) = nanmedian(ContRT(ind_Cont == 0));
%     R(2) = nanmedian(ContRT(ind_Cont == 1));
%     R(3) = nanmedian(ContRT(ind_Cont == 2));
%     R_tot = vertcat(R_tot,R);
% 
%     L(1) = nanmedian(IpsiRT(ind_Ipsi == 0));
%     L(2) = nanmedian(IpsiRT(ind_Ipsi == 1));
%     L(3) = nanmedian(IpsiRT(ind_Ipsi == 2));
%     L_tot = vertcat(L_tot,L);
%     
%     rt(end+1) = R(2)-R(1);
%     rt2(end+1) = L(2)-L(1);
%     delay(end+1) = SL(i).NormDelay;
%     
% %     % Need to split up to pre/cond/post here before removing anything
% %     ind_right = zeros(1,length(SL(i).rts_r));
% %     ind_left = zeros(1,length(SL(i).rts_l));
% %     ind_right(SL(i).righttrials(:,1)>leftbound & SL(i).righttrials(:,1)<rightbound) = 1;
% %     ind_right(SL(i).righttrials(:,1)>rightbound) = 2;
% %     ind_left(SL(i).lefttrials(:,1)>leftbound & SL(i).lefttrials(:,1)<rightbound) = 1;
% %     ind_left(SL(i).lefttrials(:,1)>rightbound) = 2;
% %     
% %     R(1) = nanmedian(SL(i).rts_r(ind_right == 0));
% %     R(2) = nanmedian(SL(i).rts_r(ind_right == 1));
% %     R(3) = nanmedian(SL(i).rts_r(ind_right == 2));
% %     R_tot = vertcat(R_tot,R);
% %     
% %     L(1) = nanmedian(SL(i).rts_l(ind_left == 0));
% %     L(2) = nanmedian(SL(i).rts_l(ind_left == 1));
% %     L(3) = nanmedian(SL(i).rts_l(ind_left == 2));
% %     L_tot = vertcat(L_tot,L);
% %     
% %     rt(end+1) = R(2)-R(1);
% %     rt2(end+1) = L(2)-L(1);
% %     delay(end+1) = SL(i).NormDelay;
%     
% end
% 
% subplot(2,2,3)
% scatter(delay,rt)
% title('Ipsi Cond - Contra Trial')
% xlabel('Normalized Stim Delay')
% ylabel('Change in Conditioning RT')
% 
% subplot(2,2,4)
% scatter(delay,rt2)
% title('Ipsi Cond - Ipsi Trial')
% xlabel('Normalized Stim Delay')
% ylabel('Change in Conditioning RT')
% 
% end