i = 10;

%% Left trials
% Aligned to trial start
Left = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); 
trig = SL(i).fs/1000*SL(i).lefttrials(:,1)';
window = [0,1]; inds = window(1)*SL(i).fs:1:window(2)*SL(i).fs;

trialinds = repmat(trig, length(inds), 1) + repmat(inds(:), 1, size(trig,1));
trialinds = floor(trialinds);
trialinds(:,trialinds(1,:)<=0) = [];
trialinds(:,trialinds(end,:)>length(SL(i).accel_raw_l)) = [];

% Raw
data = SL(i).accel_raw_l;
data = data(trialinds);
subplot(2,2,1); plot(inds/SL(i).fs,data);
title([SL(i).Animal,' Left Trial Raw, Trial Aligned'])
xlabel('Time (s)');

% Filtered 
data = u.FilterAccForMovementTimes(SL(i).accel_raw_l, SL(i).fs, 'richardson');
data = data(trialinds);
subplot(2,2,3); plot(inds/SL(i).fs,data);
title([SL(i).Animal,' Left Trial Filtered, Trial Aligned'])
xlabel('Time (s)');

% Aligned to RT
trig = SL(i).fs/1000*(SL(i).lefttrials(:,1)'+SL(i).rts_l');
trig(isnan(trig)) = [];
window = [-0.2,0.8]; inds = window(1)*SL(i).fs:1:window(2)*SL(i).fs;

trialinds = repmat(trig, length(inds), 1) + repmat(inds(:), 1, size(trig,1));
trialinds = floor(trialinds);
trialinds(:,trialinds(1,:)<=0) = [];
trialinds(:,trialinds(end,:)>length(SL(i).accel_raw_l)) = [];

% Raw
data = SL(i).accel_raw_l;
data = data(trialinds);
subplot(2,2,2); plot(inds/SL(i).fs,data);
title([SL(i).Animal,' Left Trial Raw, RT Aligned'])
xlabel('Time (s)');

% Filtered 
data = u.FilterAccForMovementTimes(SL(i).accel_raw_l, SL(i).fs, 'richardson');
data = data(trialinds);
subplot(2,2,4); plot(inds/SL(i).fs,data);
title([SL(i).Animal,' Left Trial Filtered, RT Aligned'])
xlabel('Time (s)');



%% Right trials
% Aligned to trial start
Right = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); 
Rind = [find(SL(i).righttrials(:,1) < SL(i).trig1(1),1,'last'),find(SL(i).righttrials(:,2) > SL(i).trig1(end),1)];
trig = SL(i).fs/1000*SL(i).righttrials(Rind(1):Rind(2),1)';
window = [0,1]; inds = window(1)*SL(i).fs:1:window(2)*SL(i).fs;

trialinds = repmat(trig, length(inds), 1) + repmat(inds(:), 1, size(trig,1));
trialinds = floor(trialinds);
trialinds(:,trialinds(1,:)<=0) = [];
trialinds(:,trialinds(end,:)>length(SL(i).accel_raw_l)) = [];

% Raw
data = SL(i).accel_raw_r;
data = data(trialinds);
subplot(2,2,1); plot(inds/SL(i).fs,data);
title([SL(i).Animal,' Right Trial Raw, Trial Aligned'])
xlabel('Time (s)');

% Filtered 
data = u.FilterAccForMovementTimes(SL(i).accel_raw_r, SL(i).fs, 'richardson');
data = data(trialinds);
subplot(2,2,3); plot(inds/SL(i).fs,data);
title([SL(i).Animal,' Right Trial Filtered, Trial Aligned'])
xlabel('Time (s)');

% Aligned to RT
trig = SL(i).fs/1000*(SL(i).righttrials(Rind(1):Rind(2),1)'+SL(i).rts_r(Rind(1):Rind(2))');
trig(isnan(trig)) = [];
window = [-0.2,0.8]; inds = window(1)*SL(i).fs:1:window(2)*SL(i).fs;

trialinds = repmat(trig, length(inds), 1) + repmat(inds(:), 1, size(trig,1));
trialinds = floor(trialinds);
trialinds(:,trialinds(1,:)<=0) = [];
trialinds(:,trialinds(end,:)>length(SL(i).accel_raw_l)) = [];

% Raw
data = SL(i).accel_raw_r;
data = data(trialinds);
subplot(2,2,2); plot(inds/SL(i).fs,data);
title([SL(i).Animal,' Right Trial Raw, RT Aligned'])
xlabel('Time (s)');

% Filtered 
data = u.FilterAccForMovementTimes(SL(i).accel_raw_r, SL(i).fs, 'richardson');
data = data(trialinds);
subplot(2,2,4); plot(inds/SL(i).fs,data);
title([SL(i).Animal,' Right Trial Filtered, RT Aligned'])
xlabel('Time (s)');

%% Print
file = 'C:\Users\richy.yun\Dropbox\repos\abogaard\efetz\RT manuscript\figures\Snips.ps';
print('-opengl',Left,'-dpsc2',file,'-append');
print('-opengl',Right,'-dpsc2',file,'-append');
close(Left);
close(Right);
