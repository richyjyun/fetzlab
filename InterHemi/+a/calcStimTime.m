% calcStimDelay

load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataKato.mat')

SL = SL(9);

D = SL.Date;
Sessions = strsplit(char(SL.Sessions),'_');
dwn = 1;
trig1 = []; lefttrials = []; righttrials = []; lefttrialsuccess = []; righttrialsuccess = []; data = []; last = 0; accel = [];
for i = 1:3
    trainfile = [char(D),'_',char(Sessions(i)),'.f32'];
    gugfile = [char(D),'_',char(Sessions(i))];
    if(~exist([gugfile,'.i16']) || ~exist([gugfile,'.bin']) || ~exist([gugfile,'.cfg']))
        return;
    end
    if(~exist(trainfile))
        disp(['making f32 file for ' gugfile]);
        cd('F:\Dropbox\repos\abogaard\efetz\DualAccelRt\code');
        err = utils.trainalign2(gugfile, 1);
    end
    
    [acc, t1, lt, rt, fs1, lts, rts] = u.LoadTrainKato(trainfile,dwn);
    [d, fs, chnm, ~] = u.LoadGug(gugfile, dwn);
    
    lefttrials = [lefttrials;last+lt];
    righttrials = [righttrials;last+rt];
    lefttrialsuccess = [lefttrialsuccess;lts];
    righttrialsuccess = [righttrialsuccess;rts];
    data = [data;d];
    if(i == 2)
        %             if(strcmp(SL.Condition(end),'M'))
        %                 trig1 = u.getKatoTrig(gugfile,fs,1);
        %                 trig2 = u.getKatoTrig(gugfile,fs,2);
        %             end
        if(any(SL.Condition=='I'))
            trig1 = last+rt(:,1);
            condstart = length(righttrials)-length(rt);
        elseif(any(SL.Condition=='C'))
            trig1 = last+lt(:,1);
            condstart = length(lefttrials)-length(lt);
        end
    end
    last = last+length(acc);
    accel = [accel;acc];
end
StimChns = SL.Stim_Loc;

SL.accel_raw_l = accel(:,1);
SL.accel_raw_r = accel(:,2);
SL.righttrials = righttrials;
SL.lefttrials = lefttrials;
SL.lefttrialsuccess = lefttrialsuccess;
SL.righttrialsuccess = righttrialsuccess;
SL.fs = fs;
SL = a.AppendReactionTimes(SL);

bpf = [100 500]; % chL/R transform: 1. bandpass filter (Hz) (around beta wave range)
[bbpf,abpf] = butter(1,bpf/(fs/2)); % 1st order bandpass (see daq_sapi_*)
Filter = filtfilt(bbpf,abpf,double(data));

Trials = lefttrials;
Success = lefttrialsuccess;
RT = SL.rts_l;

threshold = 50;
window = 600;
window = floor(window*SL.fs/1000);
delays = nan(1,length(trig1));
stimtime = nan(1,length(trig1));
rt = nan(1,length(trig1));
for i = 1:length(trig1)
    rt(i) = RT(i);
    if(Success(condstart+i) && ~isnan(RT(condstart+i)))
        stim = find(Filter(trig1(i):trig1(i)+window,15)>threshold,1)*1000/SL.fs;
        if(isempty(stim))
            continue;
        end
        stimtime(i) = stim;
        delays(i) = stimtime(i)-RT(i);
    end
end


window = 600;
window = floor(window*SL.fs/1000);
figure
inds = -window/3:1:window;
trialinds = repmat(trig1', length(inds), 1) + repmat(inds(:), 1, size(trig1,1));
temp = Filter(:,15);
Snips = temp(floor(trialinds));
plot(inds,Snips)

    LFilter = u.FilterAccForMovementTimes(accel(:,1), SL.fs, 'richardson');

