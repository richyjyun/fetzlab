%% Right hemisphere
% Left trial

figure; 

day = 1;
chn = 3;

subaxis(1,2,2, 'padding', .01, 'spacing', .01, 'spacingvert', .03)
win = 0.8; SL = UbiSLNeuro(day); t = SL.tneuro;
for i = 1:3
    switch i
        case 1
            c = [0,0,1];
        case 2
            c = [0.4,0.4,1];
        case 3
            c = [0.7,0.7,1];
    end
    hold on; plot(t(t>-win & t<win),squeeze(SL.Lbeta(i, chn, t>-win & t<win)),'color',c,'linewidth',2)
end
% legend('Pre','Stim','Post'); legend('boxoff'); 
axis off;

subaxis(1,2,2, 'padding', .01, 'spacing', .01, 'spacingvert', .03)
win = 0.8; SL = UbiSLNeuro(1); t = SL.tneuro;
for i = 1:3
    switch i
        case 1
            c = [1,0,0];
        case 2
            c = [1,0.4,0.4];
        case 3
            c = [1,0.7,0.7];
    end
    hold on; plot(t(t>-win & t<win),squeeze(SL.Lgamma(i, chn, t>-win & t<win)),'color',c,'linewidth',2)
end
yl1 = ylim; title('Ipsi Trials')
h(1) = plot(nan,nan,'r','color',[1,0,0],'linewidth',2);
h(2) = plot(nan,nan,'r','color',[1,0.4,0.4],'linewidth',2);
h(3) = plot(nan,nan,'r','color',[1,0.7,0.7],'linewidth',2);
legend(h,'Pre','Stim','Post'); legend('boxoff'); 
axis off;

subaxis(1,2,2, 'padding', .01, 'spacing', .01, 'spacingvert', .03)

SLInd = 2; fs= UbiSL(SLInd).fs;
trig = UbiSL(SLInd).righttrials+UbiSL(SLInd).rts_r; trig = trig(UbiSL(SLInd).righttrials(:,2)<UbiSL(SLInd).trig1(1));
Accel = u.FilterAccForMovementTimes(UbiSL(SLInd).accel_raw_r,UbiSL(SLInd).fs, 'richardson');
inds = (-win*fs:1:win*fs);
trig(isnan(trig)) = []; trig = trig*fs/1000;
trialinds = repmat(trig', length(inds), 1) + repmat(inds(:), 1, size(trig,1));
trialinds(:,trialinds(1,:)<=0) = []; trialinds(:,trialinds(end,:) > length(Accel)) = [];
Snip = Accel(floor(trialinds));
Snip = mean(Snip,2);
inds = inds/1000;

plot(inds,Snip/max(Snip)*yl1(2)/4,'k','linewidth',2);

xticks([-win,0,win]);


yl2 = ylim; ylim(yl2);


% Right trial
day = 1;
chn = 3;

subaxis(1,2,1, 'padding', .01, 'spacing', .01, 'spacingvert', .03)
win = 0.8; SL = UbiSLNeuro(day); t = SL.tneuro;
for i = 1:3
    switch i
        case 1
            c = [0,0,1];
        case 2
            c = [0.4,0.4,1];
        case 3
            c = [0.7,0.7,1];
    end
    hold on; plot(t(t>-win & t<win),squeeze(SL.Rbeta(i, chn, t>-win & t<win)),'color',c,'linewidth',2)
end
legend('Pre','Stim','Post'); legend('boxoff'); 
axis off;

subaxis(1,2,1, 'padding', .01, 'spacing', .01, 'spacingvert', .03)
win = 0.8; SL = UbiSLNeuro(1); t = SL.tneuro;
for i = 1:3
    switch i
        case 1
            c = [1,0,0];
        case 2
            c = [1,0.4,0.4];
        case 3
            c = [1,0.7,0.7];
    end
    hold on; plot(t(t>-win & t<win),squeeze(SL.Rgamma(i, chn, t>-win & t<win)),'color',c,'linewidth',2)
end
title('Contra Trials')


subaxis(1,2,1, 'padding', .01, 'spacing', .01, 'spacingvert', .03)

SLInd = 2; fs= UbiSL(SLInd).fs;
trig = UbiSL(SLInd).righttrials+UbiSL(SLInd).rts_r; trig = trig(UbiSL(SLInd).righttrials(:,2)<UbiSL(SLInd).trig1(1));
Accel = u.FilterAccForMovementTimes(UbiSL(SLInd).accel_raw_r,UbiSL(SLInd).fs, 'richardson');
inds = (-win*fs:1:win*fs);
trig(isnan(trig)) = []; trig = trig*fs/1000;
trialinds = repmat(trig', length(inds), 1) + repmat(inds(:), 1, size(trig,1));
trialinds(:,trialinds(1,:)<=0) = []; trialinds(:,trialinds(end,:) > length(Accel)) = [];
Snip = Accel(floor(trialinds));
Snip = mean(Snip,2);
inds = inds/1000;

plot(inds,Snip/max(Snip)*yl1(2)/4,'k','linewidth',2);

xticks([-win,0,win]);

ylim(yl2);

a = axes; t1 = title('Contra Hemisphere');
a.Visible = 'off'; t1.Visible = 'on';


%% Left Hemisphere
% right trials
figure;
chn = 24;
subaxis(1,2,1, 'padding', .01, 'spacing', .01, 'spacingvert', .03)
win = 0.8; SL = UbiSLNeuro(day); t = SL.tneuro;
for i = 1:3
    switch i
        case 1
            c = [0,0,1];
        case 2
            c = [0.4,0.4,1];
        case 3
            c = [0.7,0.7,1];
    end
    hold on; plot(t(t>-win & t<win),squeeze(SL.Rbeta(i, chn, t>-win & t<win)),'color',c,'linewidth',2)
end
axis off; %y1 = ylim; %mid = diff(yl)/2+yl(1); ylim([mid-diff(yl1)/2,mid+diff(yl1)/2]); 
title('Contra Trials')

subaxis(1,2,1, 'padding', .01, 'spacing', .01, 'spacingvert', .03)
win = 0.8; SL = UbiSLNeuro(1); t = SL.tneuro;
for i = 1:3
    switch i
        case 1
            c = [1,0,0];
        case 2
            c = [1,0.4,0.4];
        case 3
            c = [1,0.7,0.7];
    end
    hold on; plot(t(t>-win & t<win),squeeze(SL.Rgamma(i, chn, t>-win & t<win))*0.7+0.5,'color',c,'linewidth',2)
end
axis off; %ylim(yl2);%yl = ylim; mid = diff(yl)/2+yl(1); ylim([mid-diff(yl2)/2,mid+diff(yl2)/2]);

% h(1) = plot(nan,nan,'r','color',[1,0,0],'linewidth',2);
% h(2) = plot(nan,nan,'r','color',[1,0.4,0.4],'linewidth',2);
% h(3) = plot(nan,nan,'r','color',[1,0.7,0.7],'linewidth',2);
% legend(h,'Pre','Stim','Post'); legend('boxoff'); 
axis off;


subaxis(1,2,1, 'padding', .01, 'spacing', .01, 'spacingvert', .03)

SLInd = 2; fs= UbiSL(SLInd).fs;
trig = UbiSL(SLInd).lefttrials+UbiSL(SLInd).rts_l; trig = trig(UbiSL(SLInd).lefttrials(:,2)<UbiSL(SLInd).trig1(1));
Accel = u.FilterAccForMovementTimes(UbiSL(SLInd).accel_raw_l,UbiSL(SLInd).fs, 'richardson');
inds = (-win*fs:1:win*fs);
trig(isnan(trig)) = []; trig = trig*fs/1000;
trialinds = repmat(trig', length(inds), 1) + repmat(inds(:), 1, size(trig,1));
trialinds(:,trialinds(1,:)<=0) = []; trialinds(:,trialinds(end,:) > length(Accel)) = [];
Snip = Accel(floor(trialinds));
Snip = mean(Snip,2);
inds = inds/1000;

plot(inds,Snip/max(Snip)*yl1(2)/4,'k','linewidth',2);

xticks([-win,0,win]);

yl2 = ylim; ylim(yl2);

% left trials

chn = 24;
subaxis(1,2,2, 'padding', .01, 'spacing', .01, 'spacingvert', .03)
win = 0.8; SL = UbiSLNeuro(day); t = SL.tneuro;
for i = 1:3
    switch i
        case 1
            c = [0,0,1];
        case 2
            c = [0.4,0.4,1];
        case 3
            c = [0.7,0.7,1];
    end
    hold on; plot(t(t>-win & t<win),squeeze(SL.Lbeta(i, chn, t>-win & t<win)),'color',c,'linewidth',2)
end
axis off; %yl = ylim; mid = diff(yl)/2+yl(1); ylim([mid-diff(yl1)/2,mid+diff(yl1)/2]); 
title('Ipsi Trials')

subaxis(1,2,2, 'padding', .01, 'spacing', .01, 'spacingvert', .03)
win = 0.8; SL = UbiSLNeuro(1); t = SL.tneuro;
for i = 1:3
    switch i
        case 1
            c = [1,0,0];
        case 2
            c = [1,0.4,0.4];
        case 3
            c = [1,0.7,0.7];
    end
    hold on; plot(t(t>-win & t<win),squeeze(SL.Lgamma(i, chn, t>-win & t<win))-0.2,'color',c,'linewidth',2)
end
axis off; ylim(yl2);%yl = ylim; mid = diff(yl)/2+yl(1); ylim([mid-diff(yl2)/2,mid+diff(yl2)/2]);

% h(1) = plot(nan,nan,'r','color',[1,0,0],'linewidth',2);
% h(2) = plot(nan,nan,'r','color',[1,0.4,0.4],'linewidth',2);
% h(3) = plot(nan,nan,'r','color',[1,0.7,0.7],'linewidth',2);
% legend(h,'Pre','Stim','Post'); legend('boxoff'); axis off;


subaxis(1,2,2, 'padding', .01, 'spacing', .01, 'spacingvert', .03)

SLInd = 2; fs= UbiSL(SLInd).fs;
trig = UbiSL(SLInd).lefttrials+UbiSL(SLInd).rts_l; trig = trig(UbiSL(SLInd).lefttrials(:,2)<UbiSL(SLInd).trig1(1));
Accel = u.FilterAccForMovementTimes(UbiSL(SLInd).accel_raw_l,UbiSL(SLInd).fs, 'richardson');
inds = (-win*fs:1:win*fs);
trig(isnan(trig)) = []; trig = trig*fs/1000;
trialinds = repmat(trig', length(inds), 1) + repmat(inds(:), 1, size(trig,1));
trialinds(:,trialinds(1,:)<=0) = []; trialinds(:,trialinds(end,:) > length(Accel)) = [];
Snip = Accel(floor(trialinds));
Snip = mean(Snip,2);
inds = inds/1000;

plot(inds,Snip/max(Snip)*yl1(2)/4,'k','linewidth',2);

xticks([-win,0,win]);

a = axes; t1 = title('Ipsi Hemisphere');
a.Visible = 'off'; t1.Visible = 'on';




%% Left trial

figure; 

day = 1;
chn = 3;

subaxis(1,2,2, 'padding', .01, 'spacing', .01, 'spacingvert', .03)
win = 0.8; SL = UbiSLNeuro(day); t = SL.tneuro;
for i = 1:3
    switch i
        case 1
            c = [0,0,1];
        case 2
            c = [0.4,0.4,1];
        case 3
            c = [0.7,0.7,1];
    end
    hold on; plot(t(t>-win & t<win),squeeze(SL.Lbeta(i, chn, t>-win & t<win)),'color',c,'linewidth',2)
end
% legend('Pre','Stim','Post'); legend('boxoff'); 
axis off;

subaxis(1,2,2, 'padding', .01, 'spacing', .01, 'spacingvert', .03)
win = 0.8; SL = UbiSLNeuro(1); t = SL.tneuro;
for i = 1:3
    switch i
        case 1
            c = [1,0,0];
        case 2
            c = [1,0.4,0.4];
        case 3
            c = [1,0.7,0.7];
    end
    hold on; plot(t(t>-win & t<win),squeeze(SL.Lgamma(i, chn, t>-win & t<win)),'color',c,'linewidth',2)
end
yl1 = ylim; title('Right Hemisphere')


subaxis(1,2,2, 'padding', .01, 'spacing', .01, 'spacingvert', .03)

SLInd = 2; fs= UbiSL(SLInd).fs;
trig = UbiSL(SLInd).righttrials+UbiSL(SLInd).rts_r; trig = trig(UbiSL(SLInd).righttrials(:,2)<UbiSL(SLInd).trig1(1));
Accel = u.FilterAccForMovementTimes(UbiSL(SLInd).accel_raw_r,UbiSL(SLInd).fs, 'richardson');
inds = (-win*fs:1:win*fs);
trig(isnan(trig)) = []; trig = trig*fs/1000;
trialinds = repmat(trig', length(inds), 1) + repmat(inds(:), 1, size(trig,1));
trialinds(:,trialinds(1,:)<=0) = []; trialinds(:,trialinds(end,:) > length(Accel)) = [];
Snip = Accel(floor(trialinds));
Snip = mean(Snip,2);
inds = inds/1000;

plot(inds,Snip/max(Snip)*yl1(2)/4,'k','linewidth',2);

xticks([-win,0,win]);


yl2 = ylim; ylim(yl2);




chn = 24;
subaxis(1,2,1, 'padding', .01, 'spacing', .01, 'spacingvert', .03)
win = 0.8; SL = UbiSLNeuro(day); t = SL.tneuro;
for i = 1:3
    switch i
        case 1
            c = [0,0,1];
        case 2
            c = [0.4,0.4,1];
        case 3
            c = [0.7,0.7,1];
    end
    hold on; plot(t(t>-win & t<win),squeeze(SL.Lbeta(i, chn, t>-win & t<win)),'color',c,'linewidth',2)
end
axis off; %yl = ylim; mid = diff(yl)/2+yl(1); ylim([mid-diff(yl1)/2,mid+diff(yl1)/2]); 
title('Left Hemisphere')

subaxis(1,2,1, 'padding', .01, 'spacing', .01, 'spacingvert', .03)
win = 0.8; SL = UbiSLNeuro(1); t = SL.tneuro;
for i = 1:3
    switch i
        case 1
            c = [1,0,0];
        case 2
            c = [1,0.4,0.4];
        case 3
            c = [1,0.7,0.7];
    end
    hold on; plot(t(t>-win & t<win),squeeze(SL.Lgamma(i, chn, t>-win & t<win)),'color',c,'linewidth',2)
end
axis off; ylim(yl2);%yl = ylim; mid = diff(yl)/2+yl(1); ylim([mid-diff(yl2)/2,mid+diff(yl2)/2]);

% h(1) = plot(nan,nan,'r','color',[1,0,0],'linewidth',2);
% h(2) = plot(nan,nan,'r','color',[1,0.4,0.4],'linewidth',2);
% h(3) = plot(nan,nan,'r','color',[1,0.7,0.7],'linewidth',2);
% legend(h,'Pre','Stim','Post'); legend('boxoff'); axis off;


subaxis(1,2,1, 'padding', .01, 'spacing', .01, 'spacingvert', .03)

SLInd = 2; fs= UbiSL(SLInd).fs;
trig = UbiSL(SLInd).lefttrials+UbiSL(SLInd).rts_l; trig = trig(UbiSL(SLInd).lefttrials(:,2)<UbiSL(SLInd).trig1(1));
Accel = u.FilterAccForMovementTimes(UbiSL(SLInd).accel_raw_l,UbiSL(SLInd).fs, 'richardson');
inds = (-win*fs:1:win*fs);
trig(isnan(trig)) = []; trig = trig*fs/1000;
trialinds = repmat(trig', length(inds), 1) + repmat(inds(:), 1, size(trig,1));
trialinds(:,trialinds(1,:)<=0) = []; trialinds(:,trialinds(end,:) > length(Accel)) = [];
Snip = Accel(floor(trialinds));
Snip = mean(Snip,2);
inds = inds/1000;

plot(inds,Snip/max(Snip)*yl1(2)/4,'k','linewidth',2);

xticks([-win,0,win]);

a = axes; t1 = title('Left Trial');
a.Visible = 'off'; t1.Visible = 'on';



        
        
%% right trial

figure; 

day = 1;
chn = 3;

                subaxis(1,2,2, 'padding', .01, 'spacing', .01, 'spacingvert', .03)
win = 0.8; SL = UbiSLNeuro(day); t = SL.tneuro;
for i = 1:3
    switch i
        case 1
            c = [0,0,1];
        case 2
            c = [0.4,0.4,1];
        case 3
            c = [0.7,0.7,1];
    end
    hold on; plot(t(t>-win & t<win),squeeze(SL.Rbeta(i, chn, t>-win & t<win)),'color',c,'linewidth',2)
end
legend('Pre','Stim','Post'); legend('boxoff'); 
axis off;

                subaxis(1,2,2, 'padding', .01, 'spacing', .01, 'spacingvert', .03)
win = 0.8; SL = UbiSLNeuro(1); t = SL.tneuro;
for i = 1:3
    switch i
        case 1
            c = [1,0,0];
        case 2
            c = [1,0.4,0.4];
        case 3
            c = [1,0.7,0.7];
    end
    hold on; plot(t(t>-win & t<win),squeeze(SL.Rgamma(i, chn, t>-win & t<win)),'color',c,'linewidth',2)
end
title('Right Hemisphere')


subaxis(1,2,2, 'padding', .01, 'spacing', .01, 'spacingvert', .03)

SLInd = 2; fs= UbiSL(SLInd).fs;
trig = UbiSL(SLInd).righttrials+UbiSL(SLInd).rts_r; trig = trig(UbiSL(SLInd).righttrials(:,2)<UbiSL(SLInd).trig1(1));
Accel = u.FilterAccForMovementTimes(UbiSL(SLInd).accel_raw_r,UbiSL(SLInd).fs, 'richardson');
inds = (-win*fs:1:win*fs);
trig(isnan(trig)) = []; trig = trig*fs/1000;
trialinds = repmat(trig', length(inds), 1) + repmat(inds(:), 1, size(trig,1));
trialinds(:,trialinds(1,:)<=0) = []; trialinds(:,trialinds(end,:) > length(Accel)) = [];
Snip = Accel(floor(trialinds));
Snip = mean(Snip,2);
inds = inds/1000;

plot(inds,Snip/max(Snip)*yl1(2)/4,'k','linewidth',2);

xticks([-win,0,win]);

ylim(yl2);





chn = 24;
                subaxis(1,2,1, 'padding', .01, 'spacing', .01, 'spacingvert', .03)
win = 0.8; SL = UbiSLNeuro(day); t = SL.tneuro;
for i = 1:3
    switch i
        case 1
            c = [0,0,1];
        case 2
            c = [0.4,0.4,1];
        case 3
            c = [0.7,0.7,1];
    end
    hold on; plot(t(t>-win & t<win),squeeze(SL.Rbeta(i, chn, t>-win & t<win)),'color',c,'linewidth',2)
end
axis off; %yl = ylim; mid = diff(yl)/2+yl(1); ylim([mid-diff(yl1)/2,mid+diff(yl1)/2]); 
title('Left Hemisphere')

subaxis(1,2,1, 'padding', .01, 'spacing', .01, 'spacingvert', .03)
win = 0.8; SL = UbiSLNeuro(1); t = SL.tneuro;
for i = 1:3
    switch i
        case 1
            c = [1,0,0];
        case 2
            c = [1,0.4,0.4];
        case 3
            c = [1,0.7,0.7];
    end
    hold on; plot(t(t>-win & t<win),squeeze(SL.Rgamma(i, chn, t>-win & t<win))*0.8,'color',c,'linewidth',2)
end
axis off; ylim(yl2);%yl = ylim; mid = diff(yl)/2+yl(1); ylim([mid-diff(yl2)/2,mid+diff(yl2)/2]);

h(1) = plot(nan,nan,'r','color',[1,0,0],'linewidth',2);
h(2) = plot(nan,nan,'r','color',[1,0.4,0.4],'linewidth',2);
h(3) = plot(nan,nan,'r','color',[1,0.7,0.7],'linewidth',2);
legend(h,'Pre','Stim','Post'); legend('boxoff'); 
axis off;


subaxis(1,2,1, 'padding', .01, 'spacing', .01, 'spacingvert', .03)

SLInd = 2; fs= UbiSL(SLInd).fs;
trig = UbiSL(SLInd).lefttrials+UbiSL(SLInd).rts_l; trig = trig(UbiSL(SLInd).lefttrials(:,2)<UbiSL(SLInd).trig1(1));
Accel = u.FilterAccForMovementTimes(UbiSL(SLInd).accel_raw_l,UbiSL(SLInd).fs, 'richardson');
inds = (-win*fs:1:win*fs);
trig(isnan(trig)) = []; trig = trig*fs/1000;
trialinds = repmat(trig', length(inds), 1) + repmat(inds(:), 1, size(trig,1));
trialinds(:,trialinds(1,:)<=0) = []; trialinds(:,trialinds(end,:) > length(Accel)) = [];
Snip = Accel(floor(trialinds));
Snip = mean(Snip,2);
inds = inds/1000;

plot(inds,Snip/max(Snip)*yl1(2)/4,'k','linewidth',2);

xticks([-win,0,win]);

a = axes; t1 = title('Right Trial');
a.Visible = 'off'; t1.Visible = 'on';


        










% %% Right trial
% 
% figure; 
% 
% day = 1;
% chn = 24;
% 
% subplot(3,2,1)
% win = 0.8; SL = UbiSLNeuro(day); t = SL.tneuro;
% for i = 1:3
%     switch i
%         case 1
%             c = [0,0,1];
%         case 2
%             c = [0.4,0.4,1];
%         case 3
%             c = [0.7,0.7,1];
%     end
%     hold on; plot(t(t>-win & t<win),squeeze(SL.Rbeta(i, chn, t>-win & t<win)),'color',c,'linewidth',2)
% end
% legend('Pre','Stim','Post'); legend('boxoff'); axis off;
% yl1 = ylim;
% 
% subplot(3,2,3)
% win = 0.8; SL = UbiSLNeuro(1); t = SL.tneuro;
% for i = 1:3
%     switch i
%         case 1
%             c = [1,0,0];
%         case 2
%             c = [1,0.4,0.4];
%         case 3
%             c = [1,0.7,0.7];
%     end
%     hold on; plot(t(t>-win & t<win),squeeze(SL.Rgamma(i, chn, t>-win & t<win)),'color',c,'linewidth',2)
% end
% legend('Pre','Stim','Post'); legend('boxoff'); axis off;
% yl2 = ylim;
% 
% chn = 3;
% subplot(3,2,2)
% win = 0.8; SL = UbiSLNeuro(day); t = SL.tneuro;
% for i = 1:3
%     switch i
%         case 1
%             c = [0,0,1];
%         case 2
%             c = [0.4,0.4,1];
%         case 3
%             c = [0.7,0.7,1];
%     end
%     hold on; plot(t(t>-win & t<win),squeeze(SL.Rbeta(i, chn, t>-win & t<win)),'color',c,'linewidth',2)
% end
% axis off; yl = ylim; mid = diff(yl)/2+yl(1); ylim([mid-diff(yl1)/2,mid+diff(yl1)/2]);
% 
% subplot(3,2,4)
% win = 0.8; SL = UbiSLNeuro(1); t = SL.tneuro;
% for i = 1:3
%     switch i
%         case 1
%             c = [1,0,0];
%         case 2
%             c = [1,0.4,0.4];
%         case 3
%             c = [1,0.7,0.7];
%     end
%     hold on; plot(t(t>-win & t<win),squeeze(SL.Rgamma(i, chn, t>-win & t<win)),'color',c,'linewidth',2)
% end
% axis off; yl = ylim; mid = diff(yl)/2+yl(1); ylim([mid-diff(yl2)/2,mid+diff(yl2)/2]);
% 
% subplot(3,2,5)
