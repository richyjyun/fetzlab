% function [Sessions] = PlotControl(SL)

Cont = {};
Ipsi = {};
NTrials = [];

for i = 1:length(SL)
    cond = char(SL(i).Condition);
    
    if(~isempty(SL(i).Bad) || ~(strcmp(cond,'nostim') || strcmp(cond,'Control')...
            || (strcmp(cond,'NaN') && strcmp(SL(i).Animal,'Kato'))))
        continue;
    end
    
    if(str2num(SL(i).Date) == 20170405) %random delay stim
        continue;
    end
    
    if(strcmp(SL(i).StimHemi,'L'))
        ContRT = SL(i).rts_r;
        IpsiRT = SL(i).rts_l;
        ContT = SL(i).righttrials;
        IpsiT = SL(i).lefttrials;
    else
        ContRT = SL(i).rts_l;
        IpsiRT = SL(i).rts_r;
        ContT = SL(i).lefttrials;
        IpsiT = SL(i).righttrials;
    end

    % put trials in order
    Trials = [ContT(:,1);IpsiT(:,1)];
    RT = [ContRT;IpsiRT];
    Labels = [ones(1,length(ContT)),zeros(1,length(IpsiT))]; % 1 is contra, 0 is ipsi
    
    [Trials,order] = sort(Trials);
    RT = RT(order);
    Labels = Labels(order);
    
    %remove first 50 trials
    NTrials(end+1) = length(Trials);
    rm = 50;
    Trials = Trials(rm:end);
    RT = RT(rm:end);
    Labels = Labels(rm:end);
    
    bounds = [round(length(Trials)/3),round(2*length(Trials)/3)];
    C = find(Labels); I = find(Labels==0);
    
    Cont{1,end+1} = RT(C(C<bounds(1)));
    Cont{2,end} = RT(C(C>=bounds(1) & C< bounds(2)));
    Cont{3,end} = RT(C(C>=bounds(2)));
    
    Ipsi{1,end+1} = RT(I(I<bounds(1)));
    Ipsi{2,end} = RT(I(I>=bounds(1) & I< bounds(2)));
    Ipsi{3,end} =  RT(I(I>=bounds(2)));
end

% %% Subtracting Pre
% PreC = Cont(1,:); PreC = cellfun(@nanmedian,PreC,'UniformOutput',false);
% CondC = Cont(2,:);
% Sub = cellfun(@minus,CondC,PreC,'Un',0);
% Sub = cell2mat(Sub');
% 
% PreI = Ipsi(1,:); PreI = cellfun(@nanmedian,PreI,'UniformOutput',false);
% CondI = Ipsi(2,:);
% Sub = cellfun(@minus,CondI,PreI,'Un',0);
% Sub = cell2mat(Sub');

%% Bar graph
PreC = Cont(1,:)'; PreC = cell2mat(PreC);
CondC = Cont(2,:)'; CondC = cell2mat(CondC);
PostC = Cont(3,:)'; PostC = cell2mat(PostC);

PreI = Ipsi(1,:)'; PreI = cell2mat(PreI);
CondI = Ipsi(2,:)'; CondI = cell2mat(CondI);
PostI = Ipsi(3,:)'; PostI = cell2mat(PostI);

AllC = [PreC;CondC;PostC]; Cind = [ones(length(PreC),1);3*ones(length(CondC),1);5*ones(length(PostC),1)];
AllI = [PreI;CondI;PostI]; Iind = [2*ones(length(PreI),1);4*ones(length(CondI),1);6*ones(length(PostI),1)];
 
RT = [AllC;AllI]; idx = [Cind;Iind];

% Plot
figure; Bars = [];
nidx = length(unique(idx)); nidx = 1:(nidx/2);
positions = sort([((nidx)-1)*3+1,((nidx)-1)*3+2]);
color = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980]; % just using default matlab colors
for i = 1:length(positions)
    cind = 1;
    if mod(i,2) == 0
        cind = 2;
    end
    data = RT(idx==i);
    Bars(i) = bar(positions(i),nanmedian(data),'FaceColor',color(cind,:)); hold on;
    iqr = diff(prctile(data, [25 75]));
    err = 1.57 * iqr/sqrt(length(data));
    errorbar(positions(i),nanmedian(data),err,'k','Linewidth',1.5);
end

% legend for colors
legend( Bars(1:2),'Contra', 'Ipsi' );

% Set labels
labelpos = reshape(positions,2,length(positions)/2);
labelpos = labelpos(1,:)+diff(labelpos)/2;
set(gca,'xtick',labelpos)
set(gca,'xticklabel',labels)

% labels
title([SL(i).Animal,'\DeltaRT']); ylabel('RT (ms)'); xlabel('Stim Time (ms)')


%% Boxplot
% PreC = Cont(1,:)'; PreC = cell2mat(PreC);
% CondC = Cont(2,:)'; CondC = cell2mat(CondC);
% PostC = Cont(3,:)'; PostC = cell2mat(PostC);
% 
% PreI = Ipsi(1,:)'; PreI = cell2mat(PreI);
% CondI = Ipsi(2,:)'; CondI = cell2mat(CondI);
% PostI = Ipsi(3,:)'; PostI = cell2mat(PostI);
% 
% AllC = [PreC;CondC;PostC]; Cind = [ones(length(PreC),1);3*ones(length(CondC),1);5*ones(length(PostC),1)];
% AllI = [PreI;CondI;PostI]; Iind = [2*ones(length(PreI),1);4*ones(length(CondI),1);6*ones(length(PostI),1)];
% 
% figure; 
% RT = [AllC;AllI]; idx = [Cind;Iind];

PreC = cellfun(@nanmean,Cont(1,:));
Cont = Cont(2:3,:);
for i = 1:size(Cont,2)
    Cont{1,i} = Cont{1,i} - PreC(i);
    Cont{2,i} = Cont{2,i} - PreC(i);
end

CondC = Cont(1,:)'; CondC = cell2mat(CondC);
PostC = Cont(2,:)'; PostC = cell2mat(PostC);

PreI = cellfun(@nanmean,Ipsi(1,:));
Ipsi = Ipsi(2:3,:);
for i = 1:size(Ipsi,2)
    Ipsi{1,i} = Ipsi{1,i} - PreI(i);
    Ipsi{2,i} = Ipsi{2,i} - PreI(i);
end

CondI = Ipsi(1,:)'; CondI = cell2mat(CondI);
PostI = Ipsi(2,:)'; PostI = cell2mat(PostI);

AllC = [CondC;PostC]; Cind = [ones(length(CondC),1);3*ones(length(PostC),1)];
AllI = [CondI;PostI]; Iind = [2*ones(length(CondI),1);4*ones(length(PostI),1)];

figure; 
RT = [AllC;AllI]; idx = [Cind;Iind];

% Plot
positions = sort([1:2,(1:2)+1/3]);
boxplot(RT,idx, 'Notch','on', 'positions', positions,'Whisker',0,'Symbol','');
hold on; xl = xlim;
line(xl,[0,0],'linestyle','--','color',[0,0,0])

% Color boxes
labels = {'Cond','Post'};
color = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980]; % just using default matlab colors
color = repmat(color,length(labels),1); color = flipud(color); % need to flip because call is backwards
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),color(j,:),'FaceAlpha',.5); % set face colors
end

% legend for colors
c = get(gca, 'Children');
[~,leg] = legend(c(1:2), 'Contra', 'Ipsi' );
PatchInLegend = findobj(leg, 'type', 'patch');
set(PatchInLegend, 'facealpha', 0.5);

% labels
title([SL(i).Animal,' Control \DeltaRT']); ylabel('RT (ms)'); xlabel('Stim Time (ms)')

% Replotting to have box plot on top
boxplot(RT,idx, 'Notch','on', 'positions', positions,'Whisker',0,'Symbol','');

% Set labels
labelpos = reshape(positions,2,length(positions)/2);
labelpos = labelpos(1,:)+diff(labelpos)/2;
set(gca,'xtick',labelpos)
set(gca,'xticklabel',labels)

ylim([-150,150]);

%% Stats
str = ''; sig = 0.05;

Cidx = [1,3];
Iidx = [2,4];

Zstats = length(Cidx)+length(Iidx);
for i = 1:length(Cidx)
    Zstats(i) = signrank(RT(idx==Cidx(i)));
    n = length(RT(idx==Cidx(i)));
    stat = Zstats(i);
    sstr = 'ns';
    if(stat < sig)
        sstr = 's';
    end
    str = sprintf('%sContra%d From Zero (%d) = %0.3e\t%s\n',str,i,n,stat,sstr);
end

str = sprintf('%s\n',str);

for i = 1:length(Iidx)
    Zstats(i) = signrank(RT(idx==Iidx(i)));
    n = length(RT(idx==Iidx(i)));
    stat = Zstats(i);
    sstr = 'ns';
    if(stat < sig)
        sstr = 's';
    end
    str = sprintf('%sIpsi%d From Zero (%d) = %0.3e \t%s\n',str,i,n,stat,sstr);
end

str = sprintf('%s\n',str);

Compstats = length(Cidx);
for i = 1:length(Cidx)
    Compstats(i) = ranksum(RT(idx==Cidx(i)),RT(idx==Iidx(i)));
    stat = Compstats(i);
    sstr = 'ns';
    if(stat < sig)
        sstr = 's';
    end
    str = sprintf('%sContra%d vs Ipsi%d = %0.3e\t%s\n',str,i,i,stat,sstr);
end

str = sprintf('%s\n',str);

Cstats = nan(length(Cidx)-1);
for i = 1:length(Cidx)
    for j = i+1:length(Cidx)
        Cstats(i,j) = ranksum(RT(idx==Cidx(i)),RT(idx==Cidx(j)));
        stat = Cstats(i,j);
        sstr = 'ns';
        if(stat < sig)
            sstr = 's';
        end
        str = sprintf('%sContra%d vs Contra%d = %0.3e\t%s\n',str,i,j,stat,sstr);
    end
end

str = sprintf('%s\n',str);

Istats = nan(length(Iidx)-1);
for i = 1:length(Iidx)
    for j = i+1:length(Iidx)
        Istats(i,j) = ranksum(RT(idx==Iidx(i)),RT(idx==Iidx(j)));
        stat = Istats(i,j);
        sstr = 'ns';
        if(stat < sig)
            sstr = 's';
        end
        str = sprintf('%sIpsi%d vs Ipsi%d = %0.3e   \t%s\n',str,i,j,stat,sstr);
    end
end
str = sprintf('%s\n\n',str);


fileID = fopen('C:\Users\richy.yun\Dropbox\repos\abogaard\efetz\RT manuscript\figures\Stats_Control.txt','a');
fprintf(fileID,'%s Stats\n\n',SL(i).Animal);
fprintf(fileID,str);
fclose(fileID);










