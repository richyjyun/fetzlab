function [C] = PlotNeuralCoherence(SL)
%
% Plots coherence spectrum between several channels
%
% RJY April 12 2017

C = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]);
I = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]);

% Loading in the data. They should have the same fs.
if(strcmp(SL.Animal,'Ubi'))
    D = SL.Date;
    S = SL.Session_Guger_Train;
    Session = [char(D),'_',char(S(2))];
    trainfile = [Session,'.f32'];
    gugfile = Session;
     
    if(~exist(trainfile) || ~exist([gugfile,'.bin']) || ~exist([gugfile,'.cfg']))
        return;
    end
    
    StimChns = strsplit(SL.Stim_Loc,'/');
    n = char(StimChns(1));
    StimChns(2) = cellstr([n(1:4),char(StimChns(2))]);
    
    load([gugfile, '.cfg'],'-mat');
    
    nchans = length(find(UI.ch_enabled))+length(find(UI.ga_trigger)); % number of channels
    if isfield(UI, 'wd_channel')
        nchans = nchans + 1; % Window discriminator indicates daqdiscrim file.
    else
        nchans = nchans + 2; % Assume daqbinmanual was used.
    end
    chnm = cell(nchans,1); % name of channels
    chgu = zeros(nchans,1); % gUSBamp that channel was recorded on
    ii = 0;
    for iga = 1:length(UI.ga_trigger)
        ind = find(UI.ch_enabled(iga,:));
        chnm(ii+1:ii+length(ind)) = UI.ch_name(iga,ind)';
        chgu(ii+1:ii+length(ind)) = iga*ones(length(ind),1);
        ii = ii + length(ind);
        if UI.ga_trigger(iga), chnm{ii+1} = ['trig ' num2str(iga)]; chgu(ii+1) = iga; ii = ii + 1; end
    end
    if isfield(UI, 'wd_channel')
        chnm{end} = 'Discrim';  % One extra channel in daqdiscrim files.
    else
        chnm{end-1} = 'Behave1'; % Two extra chanels in daqbimanual files.
        chnm{end} = 'Behave2';
    end

    GoodChns = {{'RMC 3.5'},{'RMC 6.5'},{'LMC 3.5'},{'LMC 4.5'},{'LMC 8'}};
    GoodChns(strcmp([GoodChns{:}], char(StimChns(1)))|strcmp([GoodChns{:}], char(StimChns(2)))) = [];
    
    iChns = [];
    for i = 1:length(GoodChns)
        iChns(end+1) = find(strcmp(char(GoodChns{i}),chnm));
    end
    
    data = [];
    for i = 1:length(iChns)
        datfid = fopen([gugfile,'.bin'],'r');
        offset = (iChns(i)-1) * 4;
        fseek(datfid, offset, -1);
        skip_bytes = (nchans - 1) * 4;
        [samples, ~] = fread(datfid, '*single', skip_bytes);
        data(:,i) = samples;
    end
    clear samples;
    
    % down sampling rate
    dwn = 1;
    [~, trig1, ~, lefttrials, righttrials, fs, lefttrialsuccess, righttrialsuccess] = u.LoadTrainUbi(trainfile,dwn);
    
elseif(strcmp(SL.Animal,'Kato'))
    D = SL.Date;
    Sessions = strsplit(char(SL.Sessions),'_');
    dwn = 5;
    trig1 = []; lefttrials = []; righttrials = []; lefttrialsuccess = []; righttrialsuccess = []; data = []; last = 0;
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
        
        [acc, t1, lt, rt, ~, lts, rts] = u.LoadTrainKato(trainfile,dwn);
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
            elseif(any(SL.Condition=='C'))
                trig1 = last+lt(:,1);
            end
        end
        last = last+length(acc);
    end
    StimChns = SL.Stim_Loc;
elseif(strcmp(SL.Animal,'Igor'))
    D = SL.Date;
    Sessions = strsplit(char(SL.Sessions),'_');
    data = [];
    for i = 1:3
        gugfile = [char(D),'_',char(Sessions(i))];
        if( ~exist([gugfile,'.bin']) || ~exist([gugfile,'.cfg']))
            return;
        end
        [d, fs, chnm] = u.ImportIgorNeuralData(gugfile);
        
        data = [data;d];
        
    end
    StimChns = SL.Stim_Loc; trig1 = SL.trig1; 
    lefttrials = SL.lefttrials; righttrials = SL.righttrials;
    lefttrialsuccess = SL.lefttrialsuccess; righttrialsuccess = SL.righttrialsuccess;
    dwn = 1;
end

% Sometimes there are different number of trials than before due to
% downsampling. Just accounting for that.
if(~strcmp(SL.Animal,'Igor'))
    temp = lefttrials*1000/(fs*dwn);
    check = zeros(length(lefttrials),1);
    for i = 1:length(SL.lefttrials)
        norm = abs(temp(:,1)-SL.lefttrials(i,1));
        ind = find(norm == min(norm));
        check(ind) = 1;
    end
    lefttrials = lefttrials(find(check),:);
    lefttrialsuccess = lefttrialsuccess(find(check));
    
    temp = righttrials*1000/(fs*dwn);
    check = zeros(1,length(righttrials));
    for i = 1:length(SL.righttrials)
        norm = abs(temp(:,1)-SL.righttrials(i,1));
        ind = find(norm == min(norm));
        check(ind) = 1;
    end
    righttrials = righttrials(find(check),:);
    righttrialsuccess = righttrialsuccess(find(check));
end

PlotCoherence(C,lefttrials,SL.rts_l,lefttrialsuccess,trig1,fs,data);



end

function PlotCoherence(fig,trials,rt,success,trig,fs,data)
figure(fig)
window = 1.2;
window = floor(window*fs);
bounds = [trig(1),trig(end)];
rt = floor(rt./1000*fs);
for i = 1:2
    for j = 1:2
        subplot(4,1,(i-1)*2+j)
        % Grab median snippets from data(:,i) and data(:,j+2)
        pre1 = []; pre2 = []; cond1 = []; cond2=[]; post1=[]; post2=[];
        for k = 1:length(trials)
            start = trials(k,1)+rt(k)-window; fin = trials(k,1)+rt(k)+window;
            if(start<0 || fin>length(data) || isnan(rt(k)))
                continue;
            end
            if(trials(k,2) < bounds(1))
                pre1(end+1,:) = data(start:fin,i);
                pre2(end+1,:) = data(start:fin,j+2);
            elseif(trials(k,1)>(bounds(2)+50))
                post1(end+1,:) = data(start:fin,i);
                post2(end+1,:) = data(start:fin,j+2);
            elseif(~any(trig>start-50 & trig1<fin))
                cond1(end+1,:) = data(start:fin,i);
                cond2(end+1,:) = data(start:fin,j+2);
            end
        end
        pre1 = median((pre1-median(pre1,2)));
        pre2 = median((pre2-median(pre2,2)));
        cond1 = median((cond1-median(cond1,2)));
        cond2 = median((cond2-median(cond2,2)));
        post1 = median((post1-median(post1,2)));
        post2 = median((post2-median(post2,2)));
        [c1,f] = mscohere(pre1,pre2,[],[],[],fs);
        [c2,~] = mscohere(cond1,cond2,[],[],[],fs);
        [c3,~] = mscohere(post1,post2,[],[],[],fs);
        plot(f,c1,'k','LineWidth',1.5); hold on;
        plot(f,c2,'b:','LineWidth',1.5); hold on;
        plot(f,c3,'r--','LineWidth',1.5); hold off;
        
        plot(pre2,'k'); hold on;
        plot(cond2,'b:'); hold on;
        plot(post2,'r--'); hold off;
    end
end
end



