function t = PlotTwitches(SL,filter,type)

if(nargin<2)
    filter = 0;
    type = 0;
elseif(nargin<3)
    type = 0;
end

t(1) = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]);

window = 800;
row= 5;
col = 1;
pos = 1;
for i = 1:length(SL)
    if(isempty(SL(i).trig1) || ~isempty(SL(i).Bad) || strcmp(SL(i).Condition,'nostim') || ...
            strcmp(SL(i).Condition,'tonic') || strcmp(SL(i).Condition,'Control') || strcmp(SL(i).Condition,'NaN'))
        continue;
    end
    
    if(type)
        if(~strcmp(SL(i).Animal,'Ubi')) % if not ubi, assume all trials had stim in between triggers
            if strcmp(SL(i).StimHemi,'R')
                if strcmp(SL(i).Condition(1:4),'Ipsi')
                    trig = SL(i).righttrials(:,1);
                else
                    trig = SL(i).lefttrials(:,1);
                end
            else
                if strcmp(SL(i).Condition(1:4),'Ipsi')
                    trig = SL(i).lefttrials(:,1);
                else
                    trig = SL(i).righttrials(:,1);
                end
            end
            trig(trig<SL(i).trig1(1)) = [];
            trig(trig>SL(i).trig1(end)) = [];
        else % find all trials there was stimulation
            Trials = [];
            if strcmp(SL(i).Condition(1:4),'Ipsi')
                Trials = SL(i).lefttrials;
            else
                Trials = SL(i).righttrials;
            end
            stimTrial = zeros(1,length(Trials));
            for j = 1:length(SL(i).trig1)
                norm = abs(Trials(:,1) - SL(i).trig1(j));
                ind = find(norm == min(norm),1);
                stimTrial(ind) = 1;
            end
            trig = Trials(logical(stimTrial),1);
        end 
    else
        trig = SL(i).trig1;
    end
    
    inds = floor(-window*(1/4)*SL(i).fs/1000:1:window*SL(i).fs/1000);
    trialinds = repmat(trig'.*SL(i).fs./1000, length(inds), 1) + repmat(inds(:), 1, size(trig,1));
    if(~type && ~strcmp(SL(i).Animal,'Igor'))%~strcmp(SL(i).Condition(end),'M'))
        trialinds = trialinds + str2num(SL(i).Stim_Delay)*SL(i).fs/1000;
    end
    L = SL(i).accel_raw_l; R = SL(i).accel_raw_r;
    trialinds = floor(trialinds);
    trialinds(:,trialinds(1,:)<=0) = []; trialinds(:,trialinds(end,:)>min([length(L),length(R)])) = []; 
    if(filter)
        L = u.FilterAccForMovementTimes(L, SL(i).fs, 'richardson');
        R = u.FilterAccForMovementTimes(R, SL(i).fs, 'richardson');
    end
    LSnips = L(trialinds); 
    RSnips = R(trialinds); 
    if(~filter)
        LSnips = abs(LSnips - repmat(mean(LSnips,1),size(LSnips,1),1));
        RSnips = abs(RSnips - repmat(mean(RSnips,1),size(RSnips,1),1));
    end 
      
    if(pos > row*col)
        t(end+1) = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]);
        pos = 1;
    end
    subplot(row,col,pos)
    
    plot(inds./SL(i).fs.*1000, median(LSnips'),'k','LineWidth',1.5); hold on; 
    plot(inds./SL(i).fs.*1000, median(RSnips'),'r','LineWidth',1.5); hold on;
    plot([0,0],ylim,'b');
    title([SL(i).Animal,', ',num2str(SL(i).Date),', ',SL(i).StimHemi,', ',SL(i).Condition,', ',SL(i).Stim_Delay])
    set(gca,'XTick',floor(-window*(1/4)):25:floor(window)); set(gca, 'fontsize', 7);   
    xlim([floor(-window*(1/4)),floor(window)]);
    
    pos = pos+1;
end
end

