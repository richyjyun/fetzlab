%% cs_TrEN for testing transfer entropy using the TRENTOOL toolbox

close all; clear; pack

load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbi.mat')
SessionList = SL;

day = 9;

SL = SessionList(day);
disp(['Session ', SL.Date,': ',num2str(round(day/length(SessionList)*100,1)),'%']);
% if(strcmp(SL.Bad,'1') || isempty(SL.trig1) || strcmp(SL.Condition,'nostim') || strcmp(SL.Condition,'tonic') || strcmp(SL.Condition(end),'R'))
%     continue;
% end

% if(str2num(SL.Date)>=20170226 && str2num(SL.Date)<=20170306) % days with 100% stim, ignore for the time being
%     continue;
% end

%% Load Data
D = SL.Date;
S = SL.Session_Guger_Train;
Session = [char(D),'_',char(S(2))];
trainfile = [Session,'.f32'];
gugfile = Session;
% if(~exist(trainfile) || ~exist([gugfile,'.bin']) || ~exist([gugfile,'.cfg']))
%     continue;
% end

% down sampling rate
dwn = 10;
[accel, trig1, ~, lefttrials, righttrials, ~, lefttrialsuccess, righttrialsuccess] = u.LoadTrainUbi(trainfile,dwn);
[data, fs, chnm, ~] = u.LoadGug(gugfile, dwn);

if(str2num(SL.Date) == 20170127)
    trig1(1:11) = [];
end

% fo = 60;  q = 35; bw = (fo/(fs/2))/q; % comb notch filter for 60Hz noise
% [B,A] = iircomb(fs/fo,bw,'notch'); % Note type flag 'notch'
Filter = data; clear data;

% Re-calculate reaction time using upsampled, aligned data
SL.accel_raw_r = accel(:,2);
SL.accel_raw_l = accel(:,1);
SL.lefttrials = lefttrials;
SL.righttrials = righttrials;
SL.lefttrialsuccess = lefttrialsuccess;
SL.righttrialsuccess = righttrialsuccess;
SL.trig1 = trig1;
SL.fs = fs;

SL.accelfs = fs;

SL = a.AppendReactionTimes(SL);
SL = a.AppendNormalizedDelay(SL);

%% Define triggers and window
% lefttrialRT = lefttrials(:,1) + SL.rts_l./1000.*(fs*dwn);
% lefttrialsuccess(isnan(SL.rts_l)) = 0;
ofs = 250; %RP starts ~200ms before reaction time calculation

if(isempty(trig1))
    trials = [lefttrials(:,1);righttrials(:,1)];
    [~,ind] = sort(trials);
    trig1(1) = trials(ind(floor(length(ind)/4)));
    trig1(2) = trials(ind(floor(3*length(ind)/4)));
    SL.trig1 = trig1;
end

lefttrials = SL.lefttrials+(SL.rts_l-ofs)/1000*SL.fs; lefttrials = lefttrials(~isnan(SL.rts_l)& lefttrialsuccess,:);
trigPreL = lefttrials((lefttrials(:,1)<trig1(1)));
trigPostL = lefttrials((lefttrials(:,1)>trig1(end)));
trigCondL = lefttrials(lefttrials(:,1)>=trig1(1) & lefttrials(:,1)<=trig1(end),:);
remove = [];
if(strcmp(SL.Condition(1),'I'))
    for i = 1:length(trig1)
        norm = abs(trigCondL(:,1)-trig1(i));
        ind = find(norm == min(norm));
        remove(end+1) = ind;
    end
end
trigCondL(remove,:) = [];
trigCondL = trigCondL(:,1);%+SL.rts_l(length(trigPreL)+1:length(trigPreL)+length(trigCondL));
shift =0;% nanmedian(SL.rts_l);

righttrials = SL.righttrials+(SL.rts_r-ofs)/1000*SL.fs; righttrials = righttrials(~isnan(SL.rts_r) & righttrialsuccess,:);
trigPreR = righttrials(righttrials(:,1)<trig1(1),1);
trigPostR = righttrials(righttrials(:,1)>trig1(end),1);
trigCondR = righttrials(righttrials(:,1)>=trig1(1) & righttrials(:,1)<=trig1(end),:);
remove = [];
if(strcmp(SL.Condition(1),'C') && ~strcmp(SL.Condition,'Control'))
    for i = 1:length(trig1)
        norm = abs(trigCondR(:,1)-trig1(i));
        ind = find(norm == min(norm));
        remove(end+1) = ind;
    end
end
trigCondR(remove,:) = [];
trigCondR = trigCondR(:,1);

% window - time in ms after trigger
% back - time before trigger as a fraction of window
window = .5;

% Define initial channels
rchn = 6;
stim = str2num(SL.Stim_Loc(5));
if(stim == 3), lchn = 24; else, lchn = 22; end;
chnsR = [5,6];
chnsL = [lchn,lchn+1,[13,14,15,16]+17];

%% Define input variable as TRENTOOL wants it for TEprepare
% trial - {no. trials}[no. channels x no. samples] - cell array of double arrays - data for each trial
% time - {no. trials}[1 x no. samples] - cell array of double arrays (seconds) - time indeces of individual trials
% label - {no. channels x 1} - cell of strings - labels of channels
params.tapers = [3,5]; params.Fs = fs; params.fpass = [9,100]; params.trialave = 1;

trig = trigPreL;
inds = floor(-window*fs:1:0);
trialinds = floor(repmat(trig', length(inds), 1) + repmat(inds(:), 1, size(trig,1)));
trialinds(:,trialinds(1,:)<=0) = []; trialinds(:,trialinds(end,:)>length(Filter)) = [];

input = struct;
input.trial = cell(0);
for i = 1:size(trialinds,2)
    Snips = u.meanSubtract(Filter(trialinds(:,i),[chnsR,chnsL]),params);
    input.trial{i} = Snips';
end

input.time = cell(0); input.time{1} = inds/fs; input.time = repmat(input.time,1,size(trialinds,2));
input.label = chnm([chnsR,chnsL]);
input.fsample = fs;

cfg = struct; 
cfg.sgncmb = cell(0);
for i = 1:length(chnsR)
    for j = 1:length(chnsL)
        cfg.sgncmb(end+1,:) = chnm([chnsR(i),chnsL(j)])';
        cfg.sgncmb(end+1,:) = chnm([chnsL(j),chnsR(i)])';
    end
end
cfg.toi = [inds(1)/fs,inds(end)/fs];
cfg.TEcalctype = 'VW_ds';

cfg.ensemblemethod = 'no';
 
cfg.predicttime_u = 50; % u to be scanned
cfg.predicttimestepsize = 1; % time steps between u’s to be scanned

% ACT estimation and constraints on allowed ACT ( autocorelation time )
cfg.actthrvalue = 100; % threshold for ACT
cfg.maxlag = 1000;
cfg.minnrtrials = 15; % minimum acceptable number of trials

% optimizing embedding
cfg.optimizemethod ='ragwitz'; % criterion used
cfg.ragdim = 2:5; % criterion dimension
cfg.ragtaurange = [0.2 0.4]; % range for tau
cfg.ragtausteps = 5; % steps for ragwitz tau steps
cfg.repPred = 100;%size ( input . trial {1 ,1} ,2) *(3/4) ;

cfg.flagNei = 'Mass';
cfg.sizeNei = 4;

prepared = TEprepare(cfg,input);

%% Define cfg for TEsurrogatestats, and run to find TE
cfg = struct;
cfg.optdimusage = 'indivdim';
cfg.alpha = 0.05;
cfg.surrogatetype = 'trialshuffling';
cfg.extracond = 'Faes_Method';
cfg.MIcalc = 1;
cfg.shifttest = 'no';
cfg.numpermutation = 5e4;
cfg.fileidout = 'test';

results = TEsurrogatestats(cfg,prepared);

%% graph analysis
cfg = struct;
cfg.threshold = 3;
cfg.cmc = 1;
results_GA = TEgraphanalysis(cfg,results);














