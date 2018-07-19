figure;

for i = 1:length(SL)
    if(length(SL(i).Condition)<6 || ~strcmp(SL(i).Condition(1:6),'Contra') || isempty(SL(i).trig1) ||~isempty(SL(i).Bad))
        continue;
    end
    
%     stimSite = 'LMC 3/3.5';
%     if(~strcmp(SL(i).Stim_Loc,stimSite))
%         continue;
%     end
%     
    ind_var = str2num(SL(i).Stim_Amp);
    

        if(strcmp(SL(i).StimHemi,'L')) 
            contRT = SL(i).rts_r;
            ipsiRT = SL(i).rts_l;
            contT = SL(i).righttrials;
            ipsiT = SL(i).lefttrials;
        elseif(strcmp(SL(i).StimHemi,'R')) 
            contRT = SL(i).rts_l;
            ipsiRT = SL(i).rts_r;
            contT = SL(i).lefttrials;
            ipsiT = SL(i).righttrials;
        else 
            continue;
        end
 
    preI = contT(:,2)<SL(i).trig1(1);
    postI = contT(:,1)>SL(i).trig1(end);
    condI = ~(postI|preI);
    
    pre = nanmedian(contRT(preI));
    cond = nanmedian(contRT(condI));
    post = nanmedian(contRT(postI));
    
    D = char(SL(i).Date);
    subplot(2,2,1); hold on;
    scatter(ind_var,cond-pre,'k'); text(ind_var,cond-pre,D,'FontSize',7);
    subplot(2,2,3); hold on;
    scatter(ind_var,post-pre,'k'); text(ind_var,post-pre,D,'FontSize',7);
    
        
    preI = ipsiT(:,2)<SL(i).trig1(1);
    postI = ipsiT(:,1)>SL(i).trig1(end);
    condI = ~(postI|preI);
    
    pre = nanmedian(ipsiRT(preI));
    cond = nanmedian(ipsiRT(condI));
    post = nanmedian(ipsiRT(postI));
     
    subplot(2,2,2); hold on;
    scatter(ind_var,cond-pre,'k');  text(ind_var,cond-pre,D,'FontSize',7);
    subplot(2,2,4); hold on;
    scatter(ind_var,post-pre,'k'); text(ind_var,post-pre,D,'FontSize',7);
end

subplot(2,2,1); title('Cond - Pre (contra)'); xlabel('Stim Amp (uA)'); ylabel('RT (ms)')
subplot(2,2,2); title('Cond - Pre (ipsi)'); xlabel('Stim Amp (uA)'); ylabel('RT (ms)')
subplot(2,2,3); title('Post - Pre (contra)'); xlabel('Stim Amp (uA)'); ylabel('RT (ms)')
subplot(2,2,4); title('Post - Pre (ipsi)'); xlabel('Stim Amp (uA)'); ylabel('RT (ms)')

% a = axes; t1 = title(['Stim Site ',stimSite]);
% a.Visible = 'off'; t1.Visible = 'on';
    