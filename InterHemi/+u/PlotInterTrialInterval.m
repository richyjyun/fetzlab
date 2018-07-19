% load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataKatoFinal.mat');

ITI = {};
Date = {};

for i = 1:length(SL)
    if (~isempty(SL(i).Bad) || strcmp(SL(i).Condition, 'Control')...
            || strcmp(SL(i).Condition, 'NaN') || strcmp(SL(i).Condition, 'nostim')...
            || isempty(SL(i).trig1) || strcmp(SL(i).Condition(end),'R'))
        continue;
    end
    
    if(length(SL(i).Condition)<6 || ~strcmp(SL(i).Condition(1:6),'Contra'))
        continue;
    end

    Trials = [SL(i).lefttrials;SL(i).righttrials];
    RT = [SL(i).rts_l;SL(i).rts_r];
    Trials(:,2) = Trials(:,1) + RT;
    [~,idx] = sort(Trials(:,1));    
    Trials = Trials(idx,:)';
    
    idx = find(Trials(2,:) < SL(i).trig1(1),1,'last');
    Trials = Trials(:,1:idx);
    
    Trials = Trials(:);
    Trials = diff(Trials);
    Trials = Trials(2:2:length(Trials));
    
    ITI{end+1} = Trials;
    
    Date{end+1} = SL(i).Date;
    
end

all = cell2mat(ITI'); all(all>(nanmedian(all)+std(all)) | all<(nanmedian(all)-std(all))) = [];
figure; histogram(all); 
hold on; yl = ylim; plot([nanmedian(all),nanmedian(all)],yl,'r','linewidth',1.5);
xlabel('Inter Trial Interval (ms)');
ylabel('Trial Counts')
