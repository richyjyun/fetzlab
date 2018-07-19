function plotAccelData(SL,fname)

window = 0.5; % sec

h = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(h,'visible','off');

for i = 1:numel(SL)

    if(isempty(SL(i).trig1)  || ~isempty(SL(i).Bad))
        continue;
    end
    
    D = SL(i).Date;
    Session = char(D);
    disp(['Accelerometer Session ',Session])
    
    RFilter = u.FilterAccForMovementTimes(SL(i).accel_raw_r, SL(i).fs, 'richardson');
    LFilter = u.FilterAccForMovementTimes(SL(i).accel_raw_l, SL(i).fs, 'richardson');
    
    if isempty(RFilter)
        SL(i).rts_l = nan(size(SL(i).lefttrials,1),1);
        SL(i).rts_r = nan(size(SL(i).righttrials,1),1);
        continue;
    end
    
    LFilter = LFilter./SL(i).Max(1);
    RFilter = RFilter./SL(i).Max(2);
    threshold = 1/3;
    
    
    accelfs = 1000;
    if(isfield(SL(i),'accelfs'))
        accelfs = SL(i).accelfs;
    end
        
    inds = -window/2*SL(i).fs:1:window*SL(i).fs; % indeces to look backward
    stiminds = repmat(SL(i).fs/accelfs*SL(i).trig1', length(inds), 1) + repmat(inds(:), 1, size(SL(i).trig1,1));
    
    trial = []; accel = []; filt = [];
    if(length(SL(i).Condition) >= 4)
        if((strcmp(SL(i).Condition(1:4),'Ipsi') && strcmp(SL(i).StimHemi,'L')) || ...
               (strcmp(SL(i).Condition(1:6),'Contra') && strcmp(SL(i).StimHemi,'R')) )
            trial = SL(i).lefttrials;
        elseif((strcmp(SL(i).Condition(1:4),'Ipsi') && strcmp(SL(i).StimHemi,'R')) || ...
               (strcmp(SL(i).Condition(1:6),'Contra') && strcmp(SL(i).StimHemi,'L')) )
            trial = SL(i).righttrials;
        else 
            continue;
        end
    else
        continue;
    end
        
    if(strcmp(SL(i).StimHemi,'L'))
        accel = SL(i).accel_raw_r;
        filt = RFilter;
    else
        accel = SL(i).accel_raw_l;
        filt = LFilter;
    end
    
    pretrials = trial(trial(:,2)<SL(i).trig1(1),1);
    trialinds = repmat(SL(i).fs/accelfs*pretrials', length(inds), 1) + repmat(inds(:), 1, size(pretrials,1));

    stiminds = ceil(stiminds);
    trialinds = ceil(trialinds);
    trialinds(:,trialinds(1,:)<=0) = [];

%     trig = round(SL(i).fs/accelfs*SL(i).trig1);
%     stimTrials = find(trial(:,2) > trig(1) & trial(:,1)<trig(end));
%     for j = 1:length(trial)
%         if(any(trial(j,1)-10 < trig & trial(j,2)+10 > trig))
%             cond(j) = 1;
%         end
%     end
    
%     noStimCond = find(trial(:,1) < trig(1));
    
    x = inds./SL(i).fs*1000;
    subplot(2,2,1)
    plot(x,accel(trialinds)); title('Accel, Pre (trial start)'); xlim([x(1),x(end)]);
    subplot(2,2,2)
    plot(x,filt(trialinds)); title('Filtered, Pre (trial start)');xlim([x(1),x(end)]);
    
    subplot(2,2,3)
    plot(x,accel(stiminds)); title('Accel, Cond (stim trig)');xlim([x(1),x(end)]);
    subplot(2,2,4)
    plot(x,filt(stiminds)); title('Filtered, Cond (stim trig)');xlim([x(1),x(end)]);
    
    a = axes; t1 = title([SL(i).Date,' ',SL(i).Condition]);
    a.Visible = 'off'; t1.Visible = 'on';
    
    print(h, '-dpsc2', fname, '-append')
    
end

close(h);