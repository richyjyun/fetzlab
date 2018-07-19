% Plot windows of accelerometer data around the trials. Best run while
% debugging in ana.AppendReactionTimes
LFilter = u.FilterAccForMovementTimes(SL(i).accel_raw_l, SL(i).fs, 'richardson');
window = 600;
window = floor(window*SL(i).fs/1000);
figure
inds = -window/3:1:window;
trialinds = repmat(SL(i).trig1'.*SL(i).fs./1000, length(inds), 1) + repmat(inds(:), 1, size(SL(i).trig1,1));
Snips = LFilter(floor(trialinds));
plot(inds,Snips)



window = 600;
figure
accel = [];
for j = 1:length(SL(i).righttrials)%1:length(SL(22).righttrials)
    ind = SL(i).righttrials(j,2);
    if(ind<=window || ind>length(RFilter))
        continue;
    end
    hold on
    plot(RFilter(ind-window:ind))
    
    trig = (SL(i).trig1>SL(i).righttrials(j,1))&(SL(i).trig1<SL(i).righttrials(j,2));
    if(any(trig))
        ind = find(trig == 1);
        accel(end+1) = RFilter(SL(i).trig1(ind));
        pos = window - (SL(i).righttrials(j,2) - SL(i).trig1(ind));
        plot([pos,pos],[0,40])
    end
end
figure
hist(accel)


window = 600;
figure
for j = 653:length(SL(n).righttrials)
    ind = SL(n).righttrials(j,1);
    if(ind>length(RFilter)-window)
        continue;
    end
    hold on
    plot(RFilter(ind:ind+window))
end


window = 200;
figure
for j = 1:length(SL(i).lefttrials)
    ind = floor(0.001*SL(i).fs*SL(i).lefttrials(j,1));
    if(ind>length(LFilter)-window)
        continue;
    end
    hold on
    plot(LFilter(ind:ind+window))
end

% n = 1;

window = 200;
figure
accel = [];
for j = 1:length(SL(i).righttrials)
    ind = floor(0.001*SL(i).fs*SL(i).righttrials(j,2));
    if(ind<=window || ind>length(RFilter))
        continue;
    end
    hold on
    plot(RFilter(ind-window:ind))
    
        
    trig = (SL(i).triggers>SL(i).righttrials(j,2)-window/SL(i).fs*1000)&(SL(i).triggers<SL(i).righttrials(j,2));
    if(any(trig))
        ind = find(trig == 1,1);
        accel(end+1) = RFilter(floor(0.001*SL(i).fs*SL(i).triggers(ind)));
        pos = window - (floor(0.001*SL(i).fs*(SL(i).righttrials(j,2))) - floor(0.001*SL(i).fs*(SL(i).triggers(ind))));
        plot([pos,pos],[0,20])
    end
    
end
figure
hist(accel)

window = 200;
figure
for j = 1:length(SL(i).lefttrials)
    ind = floor(0.001*SL(i).fs*SL(i).lefttrials(j,2));
    if(ind<=window || ind>length(LFilter))
        continue;
    end
    hold on
    plot(LFilter(ind-window:ind))
end
nanmedian(SL(i).rts_l)