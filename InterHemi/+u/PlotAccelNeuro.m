function t = PlotAccelNeuro(SL,MRCA)

if(nargin<2 || ~strcmp(SL(1).Animal,'Ubi'))
    MRCA = 0;
end

t(1) = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]);
row = 4;
col = 1;
pos = 1;

for i = 1:length(SL)
    if(isempty(SL(i).trig1) || ~isempty(SL(i).Bad) || strcmp(SL(i).Condition,'nostim') || ...
            strcmp(SL(i).Condition,'tonic') || strcmp(SL(i).Condition,'Control') || strcmp(SL(i).Condition,'NaN'))
        continue;
    end
    
    disp(SL(i).Date);
    
    % Loading in the data. They should have the same fs.
    if(strcmp(SL(i).Animal,'Ubi'))
        D = SL(i).Date;
        S = SL(i).Session_Guger_Train;
        Session = [char(D),'_',char(S(2))];
        trainfile = [Session,'.f32'];
        gugfile = Session;
        if(~exist(trainfile) || ~exist([gugfile,'.bin']) || ~exist([gugfile,'.cfg']))
            continue;
        end
        
        StimChns = strsplit(SL(i).Stim_Loc,'/');
        n = char(StimChns(1));
        StimChns(2) = cellstr([n(1:4),char(StimChns(2))]);
        
        % down sampling rate
        dwn = 10;
        [accel, trig1, trig2, lefttrials, righttrials, ~, lefttrialsuccess, righttrialsuccess] = u.LoadTrainUbi(trainfile,dwn);
        [data, fs, chnm, ~] = u.LoadGug(gugfile, dwn); % modify so I'm only loading in one channel to make it go faster
        if(MRCA)
            chn = 12; %RMC 6.5
        else
            chn = 21; %LMC 2.5
        end
        d = data(:,chn);
        
    elseif(strcmp(SL(i).Animal,'Kato'))
        if(strcmp(SL(i).Date,'20150108')),continue;end;
        D = SL(i).Date;
        Sessions = strsplit(char(SL(i).Sessions),'_');
        dwn = 5;
        trig1 = []; lefttrials = []; righttrials = []; lefttrialsuccess = []; righttrialsuccess = []; data = []; last = 0; accel = [];
        for j = 1:3
            trainfile = [char(D),'_',char(Sessions(j)),'.f32'];
            gugfile = [char(D),'_',char(Sessions(j))];
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
            if(j == 2)
                %             if(strcmp(SL.Condition(end),'M'))
                %                 trig1 = u.getKatoTrig(gugfile,fs,1);
                %                 trig2 = u.getKatoTrig(gugfile,fs,2);
                %             end
                if(any(SL(i).Condition=='I'))
                    trig1 = last+rt(:,1);
                elseif(any(SL(i).Condition=='C'))
                    trig1 = last+lt(:,1);
                end
            end
            accel = [accel;acc];
            last = last+length(acc);
        end
        StimChns = SL(i).Stim_Loc;
        chn = 13;  % In2A seems to be giving the best artifacts
        d = data(:,chn);
        
    elseif(strcmp(SL(i).Animal,'Igor'))
        D = SL(i).Date;
        Sessions = strsplit(char(SL(i).Sessions),'_');
        data = [];
        for j = 1:3
            gugfile = [char(D),'_',char(Sessions(j))];
            if( ~exist([gugfile,'.bin']) || ~exist([gugfile,'.cfg']))
                return;
            end
            [d, fs, chnm] = u.ImportIgorNeuralData(gugfile);
            
            data = [data;d];
            
        end
        StimChns = SL(i).Stim_Loc; trig1 = SL(i).trig1;
        lefttrials = SL(i).lefttrials; righttrials = SL(i).righttrials;
        lefttrialsuccess = SL(i).lefttrialsuccess; righttrialsuccess = SL(i).righttrialsuccess;
        dwn = 1;
        StimChns = SL(i).Stim_Loc; trig1 = SL(i).trig1;
        chn = 12; % MR3A seems to be giving the best artifacts
        d = data(:,chn);
        accel = [SL(i).accel_raw_l,SL(i).accel_raw_r];
    end
    
    % i =  1;
%     clear data; 
    
    if(MRCA)
        trig = righttrials(righttrials(:,1)<trig1(1),1);
        window = 1500;
        back = 7.5;
    else
        trig = trig1;
        window = 800;
        back = 4;
    end
    
    if(str2num(D) >= 20170406)
        trig = SL(i).trig1*fs/1000;
    end
    if(strcmp(SL(i).Animal,'Igor'))
        trig = SL(i).trig1*fs/1000;
    end
    if(strcmp(SL(i).Animal,'Kato'))
        trig = trig1*fs/fs1;
    end    
    
    
    inds = floor(-window*(1/back)*fs/1000:1:window*fs/1000);
    trialinds = repmat(trig', length(inds), 1) + repmat(inds(:), 1, size(trig,1));
    if( ~strcmp(SL(i).Animal,'Igor'))%~strcmp(SL(i).Condition(end),'M'))
        trialinds = trialinds + str2num(SL(i).Stim_Delay)*fs/1000;
    end
    L = accel(:,1); R = accel(:,2);
    trialinds = floor(trialinds);
    trialinds(:,floor(trialinds(1,:)./dwn)<=0) = []; trialinds(:,trialinds(end,:)>min([length(L),length(R)])) = [];
    LSnips = L(trialinds);
    RSnips = R(trialinds);
    LSnips = abs(LSnips - repmat(mean(LSnips,1),size(LSnips,1),1));
    RSnips = abs(RSnips - repmat(mean(RSnips,1),size(RSnips,1),1));
    
    
    if(pos>row*col)
        pos = 1;
        t(end+1) = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]);
    end
    
    subplot(row,col,pos)
    plot(inds*1000/fs, median(LSnips'),'k','LineWidth',1.5); hold on;
    plot(inds*1000/fs, median(RSnips'),'r','LineWidth',1.5); hold on;
    plot([0,0],ylim,'b'); hold off;
    if(MRCA)
        RT = nanmedian(SL(i).rts_l(1:length(trig)))- str2num(SL(i).Stim_Delay);
        hold on; plot([RT,RT],ylim,'r'); hold off;
    end
    title([SL(i).Animal,', ',num2str(SL(i).Date),', ',SL(i).StimHemi,', ',SL(i).Condition,', ',SL(i).Stim_Delay])
    set(gca,'XTick',floor(-window*(1/back)):25:floor(window)); set(gca, 'fontsize', 7);
    xlim([floor(-window*(1/back)),floor(window)]);
    
    
%     inds = floor(-window*(1/back)*fs/1000:1:window*fs/1000);
%     trialinds = repmat((trig'*fs/1000), length(inds), 1) + repmat(inds(:), 1, size(trig,1));
%     if( ~strcmp(SL(i).Animal,'Igor'))%~strcmp(SL(i).Condition(end),'M'))
%         trialinds = trialinds + str2num(SL(i).Stim_Delay)*fs/1000;
%     end

    NSnips = d(floor(trialinds));
    
    
    subplot(row,col,pos+1)
    plot(inds*1000/fs,median(NSnips'-mean(NSnips',2)), 'k', 'LineWidth', 1.5); hold on;
    % plot(inds./SL(i).fs.*1000, median(NSnips'),'k','LineWidth',1.5); hold on;
%     plot(inds*1000/fs,NSnips);
    plot([0,0],ylim,'b'); hold off;
    if(MRCA)
        hold on; plot([RT,RT],ylim,'r'); hold off;
    end
    set(gca,'XTick',floor(-window*(1/back)):25:floor(window)); set(gca, 'fontsize', 7);
    xlim([floor(-window*(1/back)),floor(window)]);
    
    pos = pos+(col+1);
    
end

