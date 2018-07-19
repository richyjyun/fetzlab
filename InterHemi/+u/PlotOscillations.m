% function f = u.PlotOscillations(SL)

SLNeuro = struct([]);

for i = 1:length(SL)
    if (~isempty(SL(i).Bad) || strcmp(SL(i).Condition, 'Control')...
            || strcmp(SL(i).Condition, 'NaN') || strcmp(SL(i).Condition, 'nostim')...
            || isempty(SL(i).trig1) || strcmp(SL(i).Condition(end),'R'))
        continue;
    end
    
    D = str2num(SL(i).Date);
    Bad = 20120515;%[20120508,20120515,20120612,20120814];
    if( any(Bad==D))
        continue;
    end
    
    % load neural data
    [data, fs, chnm] = loadData(SL(i));
    
    %% Left trials
    Trials = SL(i).lefttrials; RT = SL(i).rts_l; Epochs = []; Stim = [];
    bounds = [find(Trials(:,2)<SL(i).trig1(1),1,'last'),find(Trials(:,1)>SL(i).trig1(end),1)];
    % If triggered by lefttrials for kato/igor, no conditioning epoch
    if ~strcmp(SL(i).Animal,'Ubi') && ...
            ((strcmp(SL(i).StimHemi,'L') && strcmp(SL(i).Condition(1:4),'Ipsi')) ||...
            (strcmp(SL(i).StimHemi,'R') && strcmp(SL(i).Condition(1:6),'Contra')))
        Epochs = [1,bounds(1);bounds(2),length(Trials)]; 
    else
        %%%%%%%%% NEED TO REMOVE TRIALS WHEN STIM OCCURED 
        if(strcmp(SL(i).Animal,'Ubi') && strcmp(SL(i).Condition(1:4),'Ipsi'))
            for j = 1:length(SL(i).trig1)
                norm = Trials(:,1)-((SL(i).trig1(j))+50);
                ind = find(norm>0,1)-1;
                
                if isempty(ind) || isnan(RT(ind))
                    continue;
                end

                if(abs(Trials(ind,2)-SL(i).trig1(j)) > abs(Trials(ind+1,1)-SL(i).trig1(j)))
                    ind = ind+1;
                end
                stim(end+1) = ind;
                
            end
        end
        Epochs = [1,bounds(1);bounds(1)+1,bounds(2)-1;bounds(2),length(Trials)];
    end     
    SLtemp = getOscillations(Trials,RT,Epochs,Stim,SL(i),data,fs,chnm,'L');
    
    %% Right Trials
    Trials = SL(i).righttrials; RT = SL(i).rts_r; Epochs = []; Stim = [];
    bounds = [find(Trials(:,2)<SL(i).trig1(1),1,'last'),find(Trials(:,1)>SL(i).trig1(end),1)];
    % If triggered by righttrials for kato/igor, no conditioning epoch
    if ~strcmp(SL(i).Animal,'Ubi') && ...
            ((strcmp(SL(i).StimHemi,'R') && strcmp(SL(i).Condition(1:4),'Ipsi')) ||...
            (strcmp(SL(i).StimHemi,'L') && strcmp(SL(i).Condition(1:6),'Contra')))
        Epochs = [1,bounds(1);bounds(2),length(Trials)]; 
    else
        %%%%%%%%% NEED TO REMOVE TRIALS WHEN STIM OCCURED 
        if(strcmp(SL(i).Animal,'Ubi') && strcmp(SL(i).Condition(1:4),'Contra'))
            for j = 1:length(SL(i).trig1)
                norm = Trials(:,1)-((SL(i).trig1(j))+50);
                ind = find(norm>0,1)-1;
                
                if isempty(ind) || isnan(RT(ind))
                    continue;
                end
                
                if(abs(Trials(ind,2)-SL(i).trig1(j)) > abs(Trials(ind+1,1)-SL(i).trig1(j)))
                    ind = ind+1;
                end
                stim(end+1) = ind;
                
            end
        end
        Epochs = [1,bounds(1);bounds(1)+1,bounds(2)-1;bounds(2),length(Trials)];
    end     
    
    if(isempty(SLNeuro))
        SLNeuro = getOscillations(Trials,RT,Epochs,Stim,SL(i),data,fs,chnm,'R');
    else
        SLNeuro(end+1) = getOscillations(Trials,RT,Epochs,Stim,SL(i),data,fs,chnm,'R');
    end
end

function SL = getOscillations(Trials,RT,Epochs,Stim,SL,data,fs,chnm,side)
window = 0.8; chID = []; Beta = []; Gamma = [];
for j = 1:size(Epochs,2)
    
    trial = Epochs(j,1):Epochs(j,2); trial = trial(~ismember(trial,Stim));
    trig = Trials(trial,1) + RT(trial);
    trig(isnan(trig)) = []; trig = trig*fs/1000;
    inds = (-window*fs:1:window*fs);
    trialinds = repmat(trig', length(inds), 1) + repmat(inds(:), 1, size(trig,1));
    trialinds(:,trialinds(1,:)<=0) = []; trialinds(:,trialinds(end,:) > length(data)) = [];
    
    for chn = 1:length(chnm)
        nm = chnm{chn}; 
        if(strcmp(SL.Animal,'Ubi') && ~(strcmp(nm(1),'R') || strcmp(nm(1),'L')))
             continue;
        end
       
        disp(nm);
        
        if(j==1)
            chID(end+1) = chn;
        end
        
        d = data(:,chn);
        params.tapers = [3,5]; params.Fs = fs; params.fpass = [9,100]; params.trialave = 1;
        NSnips = u.meanSubtract(d(floor(trialinds)),params);
        NSnips(:,any(NSnips > mean(std(NSnips))*10)) = []; %remove ones with outliers (noise, artifact, etc)
        
        beta = [15,30]; gamma = [50,100];
        
        %         [t,S] = getPower(beta,params);
        movingwin = [0.25,0.01]; params.fpass = beta;
        [S,t,~] = mtspecgramc(NSnips,movingwin,params);
        S = mean(S');
        Beta(j,chn,:) = S;
        
        %         [t,S] = getPower(gamma,params);
        movingwin = [0.25,0.01]; params.fpass = gamma;
        [S,t,~] = mtspecgramc(NSnips,movingwin,params);
        S = mean(S');
        Gamma(j,chn,:) = S;

    end
    
end

Beta(:,Beta(1,:,1) == 0,:) = [];
Gamma(:,Gamma(1,:,1) == 0,:) = [];

SL.tneuro = t;
SL.([side,'beta']) = Beta;
SL.([side,'gamma']) = Gamma;
SL.chnm = chnm;
SL.chID = chID;

end


function [t,S] = getPower(bp,params)
% with chronux
movingwin = [0.25,0.01]; params.fpass = bp;
[S,t,~] = mtspecgramc(NSnips,movingwin,params);
S = mean(S');
end

function [data, fs, chnm] = loadData(SL)
% Loading in the neural data. No need to load in .f32 since everything
    % should already be set to the aligned data
    if(strcmp(SL.Animal,'Ubi'))
        D = SL.Date;
        S = SL.Session_Guger_Train;
        Session = [char(D),'_',char(S(2))];
        gugfile = Session;
        if(~exist([gugfile,'.bin']) || ~exist([gugfile,'.cfg']))
            continue;
        end
        % down sampling rate
        dwn = 10;
        [data, fs, chnm, ~] = u.LoadGug(gugfile, dwn); % modify so I'm only loading in one channel to make it go faster
        
    elseif(strcmp(SL.Animal,'Kato'))
        if(strcmp(SL.Date,'20150108')),continue;end;
        D = SL.Date;
        Sessions = strsplit(char(SL.Sessions),'_');
        dwn = 5;
        data = [];
        for j = 1:3
            gugfile = [char(D),'_',char(Sessions(j))];
            if(~exist([gugfile,'.bin']) || ~exist([gugfile,'.cfg']))
                return;
            end
            [d, fs, chnm, ~] = u.LoadGug(gugfile, dwn);
            data = [data;d]; 
        end
        
    elseif(strcmp(SL.Animal,'Igor'))
        D = SL.Date;
        Sessions = strsplit(char(SL.Sessions),'_');
        data = [];
        for j = 1:3
            gugfile = [char(D),'_',char(Sessions(j))];
            if( ~exist([gugfile,'.bin']) || ~exist([gugfile,'.cfg']))
                return;
            end
            [d, fs, chnm] = u.ImportIgorNeuralData(gugfile);
            
            data = [data;d];
        end

    end
end

% 
% function f = PlotTrial(Trials,RT,Epochs,SL,data,fs,chnm,accel)
% window = 0.8;
% f(1) = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(f(1),'visible','off');
% for j = 1:size(Epochs,2)
%     switch j
%         case 1 
%             c = 0.7;
%         case 2
%             c = 0.4;
%         case 3 
%             c = 0.1;
%     end
%     
%     trig = Trials(Epochs(j,1):Epochs(j,2),1) + RT(Epochs(j,1):Epochs(j,2));
%     trig(isnan(trig)) = []; trig = trig*fs/1000;
%     inds = (-window*fs:1:window*fs);
%     trialinds = repmat(trig', length(inds), 1) + repmat(inds(:), 1, size(trig,1));
%     trialinds(:,trialinds(1,:)<=0) = []; trialinds(:,trialinds(end,:) > length(data)) = [];
%     
%     trig = Trials(Epochs(j,1):Epochs(j,2),1) + RT(Epochs(j,1):Epochs(j,2));
%     trig(isnan(trig)) = []; trig = trig*SL.fs/1000;
%     ainds = (-window*SL.fs:1:window*SL.fs);
%     accelinds = repmat(trig', length(inds), 1) + repmat(inds(:), 1, size(trig,1));
%     accelinds(:,accelinds(1,:)<=0) = []; accelinds(:,accelinds(end,:) > length(accel)) = [];
%     
%     figInd = 1;
%     for chn = 1:length(chnm)
%         nm = chnm{chn}; trial = [];
%         if(strcmp(SL.Animal,'Ubi'))
%             if(strcmp(nm(1),'R'))
%                 trial = 'Right Trials';
%             elseif(strcmp(nm(1),'L'))
%                 trial = 'Left Trials';
%             else
%                 continue;
%             end
%         end
%         disp(nm);
%         
%         d = data(:,chn);
%         params.tapers = [3,5]; params.Fs = fs; params.fpass = [9,100]; params.trialave = 1;
%         NSnips = u.meanSubtract(d(floor(trialinds)),params);
%         NSnips(:,any(NSnips > mean(std(NSnips))*10)) = []; %remove ones with outliers (noise, artifact, etc)
%         
%         if(figInd > 36)
%             f(end+1) = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(f(end),'visible','off')
%             figInd = 1;
%         end
%         
%         subaxis(6, 6, figInd, 'padding', .01, 'spacing', .01, 'spacingvert', .03)
%         
%         beta = [15,30]; gamma = [50,100];
%         
%         %         [t,S] = getPower(beta,params);
%         movingwin = [0.25,0.01]; params.fpass = beta;
%         [S,t,~] = mtspecgramc(NSnips,movingwin,params);
%         S = mean(S');
%         yyaxis left; hold on; plot(t-window,S,'color',[c,c,1],'Linewidth',2);
%         Acc = mean(accel(floor(accelinds)),2);
%         if(j==size(Epochs,2))
%             yl = ylim; hold on; plot(ainds/1000,Acc/max(Acc)*yl(2)/5,'k','Linewidth',1.5)
%         end
%         set(gca, 'fontsize', 7)
% 
%         
%         %         [t,S] = getPower(gamma,params);
%         movingwin = [0.25,0.01]; params.fpass = gamma;
%         [S,t,~] = mtspecgramc(NSnips,movingwin,params);
%         S = mean(S');
%         yyaxis right; hold on; plot(t-window,S,'color',[1,c,c],'Linewidth',2); 
%         title([chnm(chn),' ',trial],'fontsize',7)
%         xlim([t(1)-window,t(end)-window]); set(gca, 'fontsize', 7)
%         
%         %%%%%%%%%%%%  MAKE SURE TO SAVE CURVES TO STRUCT
%         %%%%%%%%%%%%  (Lbeta,Lgamma,Rbeta,Rgamma,tneuro,chnm,chind)
% %         movingwin = [0.25,0.01];
% %         [S,t,f] = mtspecgramc(NSnips,movingwin,params);
% %         figure;
% %         imagesc(t,f,log(S'))
% %         set(gca,'YDir','normal')
% %         
% %         Wn_theta = [bp(1)/(fs/2) bp(2)/(fs/2)]; % normalized by the nyquist frequency
% %         [btheta,atheta] = butter(3,Wn_theta);
% %         sig = filtfilt(btheta,atheta,NSnips');
% %         
% %         avp = mean(sig.^2);
% %         sig = abs(hilbert(sig));
% %         hold on; plot(smooth(mean(sig),20))
%         
%         figInd = figInd+1;
%     end
%     
% end
% end






