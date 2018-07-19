function h = ShowMovementTimesStats(SL)
%
% This function will generate a distribution of movement times for each
% monkey
%
% arb 15 april 2015

% Modified for Ubi, 14 March 2017

tf = cellfun(@isempty, {SL(:).rts_l}, 'unifo', 1); % there are no reaction times

SL(tf) = []; % delete those

monkeys = {'Ubi','Kato','Igor'};

slmonkeys = {SL(:).Animal};

SL = u.SortSL(SL,'Date');

h = figure;%('PaperPosition', [.2 .2 8.3 10.8], 'position', [0 0 850 1100]);

bins = 0:5:1000;

for i = 1:length(monkeys)
    
    tf = strcmp(slmonkeys, monkeys{i}); % tf where SL is the monkey of interest
    if ~any(tf), continue; end
    where = find(tf);
    
    dist_left = cellfun(@(x, y) histc(x, y), {SL(tf).rts_l}, repmat({bins}, 1, sum(tf)), 'unif', 0);
    dist_right = cellfun(@(x, y) histc(x, y), {SL(tf).rts_r}, repmat({bins}, 1, sum(tf)), 'unif', 0);
    
    dist_left = horzcat(dist_left{:});
    dist_right = horzcat(dist_right{:});
   
    subaxis(1,2,2, 'padding', .01, 'spacing', .01, 'spacingvert', .03);    imagesc(1:size(dist_left,2), bins, dist_left, [0 .6*max(dist_left(:))]), title([monkeys{i} ' Left Reaction Times'])
    set(gca,'YDir','normal'); c1 = caxis; c = colorbar; c.Label.String = 'Number of Trials';
    subaxis(1,2,1, 'padding', .01, 'spacing', .01, 'spacingvert', .03);
    imagesc(1:size(dist_right,2), bins, dist_right, [0 .6*max(dist_right(:))]), title([monkeys{i} ' Right Reaction Times'])
    ylabel('time (ms)'); set(gca,'YDir','normal'); c2 = caxis; 
    
    subaxis(1,2,2, 'padding', .01, 'spacing', .01, 'spacingvert', .03);
    caxis([min(c1(1),c2(1)),max(c1(2),c2(2))]);
    u.xticklabel_rotate(round(linspace(1,length(where),5)),30, cellfun(@num2str, {SL(where(round(linspace(1,length(where),5)))).Date}, 'unif', 0), 'interpreter','none','fontsize',7);
    subaxis(1,2,1, 'padding', .01, 'spacing', .01, 'spacingvert', .03); caxis([min(c1(1),c2(1)),max(c1(2),c2(2))]);
    u.xticklabel_rotate(round(linspace(1,length(where),5)),30, cellfun(@num2str, {SL(where(round(linspace(1,length(where),5)))).Date}, 'unif', 0), 'interpreter','none','fontsize',7);

end

print(h, '-dpsc2', 'test.ps')

%% this stuff calculates the movement detected relative to trial end
%     leftdetect = -(diff(vertcat(SL(tf).lefttrials), 1, 2) - vertcat(SL(tf).rts_l));
%     rightdetect = -(diff(vertcat(SL(tf).righttrials), 1, 2) - vertcat(SL(tf).rts_r));
%      
%     leftdetect = sort(leftdetect, 'ascend');
%     rightdetect = sort(rightdetect, 'ascend');
%     
%     leftdetect(isnan(leftdetect)) = [];
%     rightdetect(isnan(rightdetect)) = [];
%     
%     leftthresh = leftdetect(round(0.025*length(leftdetect)));
%     rightthresh = rightdetect(round(0.025*length(rightdetect)));
%     
%     % throw out bottom 2.5%?
%     figure
%     hist(leftdetect, 1000), hold on;
%     plot([leftthresh leftthresh], [0 2], '-r')
%     figure
%     hist(rightdetect, 1000), hold on;
%     plot([rightthresh rightthresh], [0 2], '-r')
%     keyboard