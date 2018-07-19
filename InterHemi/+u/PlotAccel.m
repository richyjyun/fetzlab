fname = fullfile('F:\Dropbox\repos\abogaard\efetz\U\packets', 'KatoAccel.ps');
% Delete previous ps file so it doesn't keep appending
if(exist(fname))
    delete(fname);
end

for i = 1:length(SL)
    if(strcmp(SL(i).Bad,'1') || isempty(SL(i).trig1) || strcmp(SL(i).Condition,'Control') || strcmp(SL(i).Condition,'tonic') || strcmp(SL(i).Condition(end),'R'))
        continue;
    end
    
    D = SL(i).Date;
    Session = char(D);
    disp(['Accel Session ',Session])
    
    bound(1) = find(SL(i).lefttrials(:,2)<SL(i).trig1(1),1,'last');
    bound(2) = find(SL(i).lefttrials(:,1)>SL(i).trig1(end),1);
    
    f = plotTrials(SL(i).lefttrials,SL(i).rts_l,SL(i).accel_raw_l,bound,SL(i).fs);
    a = axes; t1 = title([SL(i).Date,' ',SL(i).Condition,' Left']);
    a.Visible = 'off'; t1.Visible = 'on';
    print(f, '-dpsc2', fname, '-append')
    close(f)
    
    f = plotTrials(SL(i).righttrials,SL(i).rts_r,SL(i).accel_raw_r,bound,SL(i).fs);
    a = axes; t1 = title([SL(i).Date,' ',SL(i).Condition,' Right']);
    a.Visible = 'off'; t1.Visible = 'on';
    print(f, '-dpsc2', fname, '-append')
    close(f)
    
end


function f = plotTrials(Trials,RT,Accel,bounds,fs)
f = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(f,'visible','off');
Filter = u.FilterAccForMovementTimes(Accel, fs, 'richardson');

bounds = [1,bounds(1);bounds(1)+1,bounds(2)-1;bounds(2),length(Trials)];
for i = 1:3
    trig = Trials(bounds(i,1):bounds(i,2),1) + RT(bounds(i,1):bounds(i,2));
    trig(isnan(trig)) = []; trig = trig*fs/1000;
    window = 350;
    inds = floor(-window*fs/1000:1:window*fs/1000);
    trialinds = repmat(trig', length(inds), 1) + repmat(inds(:), 1, size(trig,1));
    
    x = inds*1000/fs;
    
    subplot(3,2,i*2-1)
    plot(x,Accel(floor(trialinds))); 
    hold on; plot(x,mean(Accel(floor(trialinds))'),'k','Linewidth',2); xlim([x(1),x(end)])

    subplot(3,2,i*2)
    plot(x,Filter(floor(trialinds))); 
    hold on; plot(x,mean(Filter(floor(trialinds))'),'k','Linewidth',2); xlim([x(1),x(end)])
end

end