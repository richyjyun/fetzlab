function  [ipsiCSipsiRT, ipsiCScontRT, contCSipsiRT, contCScontRT, ipsiCScontRTLFPbeta, ipsiCSipsiRTLFPbeta, contCScontRTLFPbeta, contCSipsiRTLFPbeta, controlRTL, ipsiCScontbeta, ipsiCSipsibeta, contCScontbeta, contCSipsibeta, controlRTR, ipsiCScontRT_rtvt, ipsiCSipsiRT_rtvt, contCScontRT_rtvt, contCSipsiRT_rtvt] = GetStats(SL, ch_hem)
%
% ch_hem is either L or R for where lfp channel is
%
% This function will gather stats for four conditions:
% 1) IPSI CS and IPSI RT
% 2) IPSI CS and CONT RT
% 3) CONT CS and IPSI RT
% 4) CONT CS and CONT RT
%

tf_LHEM = ch_hem=='L';
tf_RHEM = ch_hem=='R';

ipsiCSipsiRT = [];
ipsiCScontRT = [];
contCSipsiRT = [];
contCScontRT = [];

controlRTL = [];
controlRTR = [];

controlbeta = [];
ipsiCScontbeta = [];
ipsiCSipsibeta = [];
contCScontbeta = [];
contCSipsibeta = [];

ipsiCScontRTLFPbeta = [];
ipsiCSipsiRTLFPbeta = [];
contCScontRTLFPbeta = [];
contCSipsiRTLFPbeta = [];

ipsiCScontRTLFPgamma = [];
ipsiCSipsiRTLFPgamma = [];
contCScontRTLFPgamma = [];
contCSipsiRTLFPgamma = [];

ipsiCScontRT_rtvt = [];
ipsiCSipsiRT_rtvt = [];
contCScontRT_rtvt = [];
contCSipsiRT_rtvt = [];

for i = 1:size(SL,1) % go through all experiment sets
    
    if any(cellfun(@isempty, {SL(i,:).rts_l}, 'unif', 1)), disp(['Skipping due to empty data: ' num2str(SL(i,1).expid)]); continue; end
    
    if strcmp(SL(i,1).monkey, 'Kato') % all Kato stim hem was R
        stimhem = 'R';
    else
        if any(SL(i,1).stimchannels{1}=='L'), stimhem = 'L';
        elseif any(SL(i,1).stimchannels{1}=='R'), stimhem = 'R';
        else stimhem = 'X';
        end
    end
    
    if any(SL(i,1).cond=='L'), cond = 'L';
    elseif any(SL(i,1).cond=='R'), cond = 'R';
    else
        if strcmp(SL(i,1).cond, 'nostim')
            disp('found one')
            controlRTL(end+1,1:3) = cellfun(@(c) nanmedian(c), {SL(i, :).rts_l}, 'unif', 1);
            controlRTR(end+1,1:3) = cellfun(@(c) nanmedian(c), {SL(i, :).rts_r}, 'unif', 1);
            controlbeta = vertcat(controlbeta, vertcat(SL(i,:).avbeta)');
            
            continue
            
        end
            
        disp(['Skipped condition named ' SL(i,1).cond]); continue %error('somethin wrong 2');
    end

    switch [stimhem cond]
        
        case 'LL' % ipsilateral left stim
            ipsiCSipsiRT(end+1,1:3) = cellfun(@(c) nanmedian(c), {SL(i, :).rts_l}, 'unif', 1);
            ipsiCSipsiRT(end,4) = SL(i,2).lag_stimtarget; % append the stim lag
            ipsiCScontRT(end+1,1:3) = cellfun(@(c) nanmedian(c), {SL(i, :).rts_r}, 'unif', 1);
            ipsiCScontRT(end,4) = SL(i,2).lag_stimtarget; % append the stim lag
            
            ipsiCScontbeta = vertcat(ipsiCScontbeta, horzcat(SL(i,1).avbeta(tf_RHEM)', SL(i,2).avbeta(tf_RHEM)', SL(i,3).avbeta(tf_RHEM)'));
            ipsiCSipsibeta = vertcat(ipsiCSipsibeta, horzcat(SL(i,1).avbeta(tf_LHEM)', SL(i,2).avbeta(tf_LHEM)', SL(i,3).avbeta(tf_LHEM)'));
            
            ipsiCSipsiRTLFPbeta(:,1:6,end+1) = [mean(SL(i,1).TA_Lbeta(:,tf_LHEM),2), mean(SL(i,1).TA_Lbeta(:,tf_RHEM),2),... [pre ipsi pre cont con ipsi con cont..
                                                mean(SL(i,2).TA_Lbeta(:,tf_LHEM),2), mean(SL(i,2).TA_Lbeta(:,tf_RHEM),2),...
                                                mean(SL(i,3).TA_Lbeta(:,tf_LHEM),2), mean(SL(i,3).TA_Lbeta(:,tf_RHEM),2)];
            ipsiCScontRTLFPbeta(:,1:6,end+1) = [mean(SL(i,1).TA_Rbeta(:,tf_LHEM),2), mean(SL(i,1).TA_Rbeta(:,tf_RHEM),2),...
                                                mean(SL(i,2).TA_Rbeta(:,tf_LHEM),2), mean(SL(i,2).TA_Rbeta(:,tf_RHEM),2),...
                                                mean(SL(i,3).TA_Rbeta(:,tf_LHEM),2), mean(SL(i,3).TA_Rbeta(:,tf_RHEM),2)];    
                                            
            ipsiCScontRT_rtvt = vertcat(ipsiCScontRT_rtvt, [SL(i,2).righttrials(:,1)+SL(i,2).rts_r(:),SL(i,2).rts_r(:)]);
            ipsiCSipsiRT_rtvt = vertcat(ipsiCSipsiRT_rtvt, [SL(i,2).lefttrials(:,1)+SL(i,2).rts_l(:),SL(i,2).rts_l(:)]);
            
        case 'LR' % contralateral left stim
            contCSipsiRT(end+1,1:3) = cellfun(@(c) nanmedian(c), {SL(i, :).rts_l}, 'unif', 1);
            contCSipsiRT(end,4) = SL(i,2).lag_stimtarget; % append the stim lag
            contCScontRT(end+1,1:3) = cellfun(@(c) nanmedian(c), {SL(i, :).rts_r}, 'unif', 1);
            contCScontRT(end,4) = SL(i,2).lag_stimtarget; % append the stim lag
            
            contCScontbeta = vertcat(contCScontbeta, horzcat(SL(i,1).avbeta(tf_RHEM)', SL(i,2).avbeta(tf_RHEM)', SL(i,3).avbeta(tf_RHEM)'));
            contCSipsibeta = vertcat(contCSipsibeta, horzcat(SL(i,1).avbeta(tf_LHEM)', SL(i,2).avbeta(tf_LHEM)', SL(i,3).avbeta(tf_LHEM)'));
            
            contCSipsiRTLFPbeta(:,1:6,end+1) = [mean(SL(i,1).TA_Lbeta(:,tf_LHEM),2), mean(SL(i,1).TA_Lbeta(:,tf_RHEM),2),...
                                                mean(SL(i,2).TA_Lbeta(:,tf_LHEM),2), mean(SL(i,2).TA_Lbeta(:,tf_RHEM),2),...
                                                mean(SL(i,3).TA_Lbeta(:,tf_LHEM),2), mean(SL(i,3).TA_Lbeta(:,tf_RHEM),2)];
            contCScontRTLFPbeta(:,1:6,end+1) = [mean(SL(i,1).TA_Rbeta(:,tf_LHEM),2), mean(SL(i,1).TA_Rbeta(:,tf_RHEM),2),...
                                                mean(SL(i,2).TA_Rbeta(:,tf_LHEM),2), mean(SL(i,2).TA_Rbeta(:,tf_RHEM),2),...
                                                mean(SL(i,3).TA_Rbeta(:,tf_LHEM),2), mean(SL(i,3).TA_Rbeta(:,tf_RHEM),2)];
                                            
            contCScontRT_rtvt = vertcat(contCScontRT_rtvt, [SL(i,2).righttrials(:,1)+SL(i,2).rts_r(:),SL(i,2).rts_r(:)]);
            contCSipsiRT_rtvt = vertcat(contCSipsiRT_rtvt, [SL(i,2).lefttrials(:,1)+SL(i,2).rts_l(:),SL(i,2).rts_l(:)]);
            
        case 'RR' % ipsilateral right stim
            ipsiCSipsiRT(end+1,1:3) = cellfun(@(c) nanmedian(c), {SL(i, :).rts_r}, 'unif', 1);
            ipsiCSipsiRT(end,4) = SL(i,2).lag_stimtarget; % append the stim lag
            ipsiCScontRT(end+1,1:3) = cellfun(@(c) nanmedian(c), {SL(i, :).rts_l}, 'unif', 1);
            ipsiCScontRT(end,4) = SL(i,2).lag_stimtarget; % append the stim lag
            
            ipsiCScontbeta = vertcat(ipsiCScontbeta, horzcat(SL(i,1).avbeta(tf_LHEM)', SL(i,2).avbeta(tf_LHEM)', SL(i,3).avbeta(tf_LHEM)'));
            ipsiCSipsibeta = vertcat(ipsiCSipsibeta, horzcat(SL(i,1).avbeta(tf_RHEM)', SL(i,2).avbeta(tf_RHEM)', SL(i,3).avbeta(tf_RHEM)'));
            
            ipsiCSipsiRTLFPbeta(:,1:6,end+1) = [mean(SL(i,1).TA_Rbeta(:,tf_RHEM),2), mean(SL(i,1).TA_Rbeta(:,tf_LHEM),2),...
                                                mean(SL(i,2).TA_Rbeta(:,tf_RHEM),2), mean(SL(i,2).TA_Rbeta(:,tf_LHEM),2),...
                                                mean(SL(i,3).TA_Rbeta(:,tf_RHEM),2), mean(SL(i,3).TA_Rbeta(:,tf_LHEM),2)];
            ipsiCScontRTLFPbeta(:,1:6,end+1) = [mean(SL(i,1).TA_Lbeta(:,tf_RHEM),2), mean(SL(i,1).TA_Lbeta(:,tf_LHEM),2),...
                                                mean(SL(i,2).TA_Lbeta(:,tf_RHEM),2), mean(SL(i,2).TA_Lbeta(:,tf_LHEM),2),...
                                                mean(SL(i,3).TA_Lbeta(:,tf_RHEM),2), mean(SL(i,3).TA_Lbeta(:,tf_LHEM),2)]; 
                                            
            contCSipsiRT_rtvt = vertcat(contCSipsiRT_rtvt, [SL(i,2).righttrials(:,1)+SL(i,2).rts_r(:),SL(i,2).rts_r(:)]);
            contCScontRT_rtvt = vertcat(contCScontRT_rtvt, [SL(i,2).lefttrials(:,1)+SL(i,2).rts_l(:),SL(i,2).rts_l(:)]);                                
            
        case 'RL' % contralateral right stim
            contCSipsiRT(end+1,1:3) = cellfun(@(c) nanmedian(c), {SL(i, :).rts_r}, 'unif', 1);
            contCSipsiRT(end,4) = SL(i,2).lag_stimtarget; % append the stim lag
            contCScontRT(end+1,1:3) = cellfun(@(c) nanmedian(c), {SL(i, :).rts_l}, 'unif', 1);
            contCScontRT(end,4) = SL(i,2).lag_stimtarget; % append the stim lag
            
            contCScontbeta = vertcat(contCScontbeta, horzcat(SL(i,1).avbeta(tf_LHEM)', SL(i,2).avbeta(tf_LHEM)', SL(i,3).avbeta(tf_LHEM)'));
            contCSipsibeta = vertcat(contCSipsibeta, horzcat(SL(i,1).avbeta(tf_RHEM)', SL(i,2).avbeta(tf_RHEM)', SL(i,3).avbeta(tf_RHEM)'));
            
            contCSipsiRTLFPbeta(:,1:6,end+1) = [mean(SL(i,1).TA_Rbeta(:,tf_RHEM),2), mean(SL(i,1).TA_Rbeta(:,tf_LHEM),2),...
                                                mean(SL(i,2).TA_Rbeta(:,tf_RHEM),2), mean(SL(i,2).TA_Rbeta(:,tf_LHEM),2),...
                                                mean(SL(i,3).TA_Rbeta(:,tf_RHEM),2), mean(SL(i,3).TA_Rbeta(:,tf_LHEM),2)];
            contCScontRTLFPbeta(:,1:6,end+1) = [mean(SL(i,1).TA_Lbeta(:,tf_RHEM),2), mean(SL(i,1).TA_Lbeta(:,tf_LHEM),2),...
                                                mean(SL(i,2).TA_Lbeta(:,tf_RHEM),2), mean(SL(i,2).TA_Lbeta(:,tf_LHEM),2),...
                                                mean(SL(i,3).TA_Lbeta(:,tf_RHEM),2), mean(SL(i,3).TA_Lbeta(:,tf_LHEM),2)];
                                            
            contCScontRT_rtvt = vertcat(contCScontRT_rtvt, [SL(i,2).righttrials(:,1)+SL(i,2).rts_r(:),SL(i,2).rts_r(:)]);
            contCSipsiRT_rtvt = vertcat(contCSipsiRT_rtvt, [SL(i,2).lefttrials(:,1)+SL(i,2).rts_l(:),SL(i,2).rts_l(:)]);                                
            
    end
end