% function PlotBetaRebound(SLNeuro,SL,fname)

% fname = fullfile('F:\Dropbox\repos\abogaard\efetz\U\packets', 'BetaGammaUbi.ps');
% if(exist(fname))
%     delete(fname);
% end

dates = extractfield(SL,'Date'); window = 0.9;
Beta = struct('C',[],'I',[]); 
StimHemi = []; 
for i = 1:length(SLNeuro)
    inds = SLNeuro(i).tneuro > 0 & SLNeuro(i).tneuro < window;
    
    %% find proper SL, get change in RT
    % Doing post - pre for igor, can change to cond-pre for ubi
    SLInd = find(strcmp(dates,SLNeuro(i).Date)); trig = SL(SLInd).trig1;
    D = str2num(SL(SLInd).Date);
    if(~strcmp(SL(SLInd).Condition(1:6),'Contra'))
        continue;
    end    
 
    if(strcmp(SL(SLInd).StimHemi,'L'))
        Beta.C(end+1,:,:) = max(SLNeuro(i).Rbeta(:,:,inds),[],3);
        Beta.I(end+1,:,:) = max(SLNeuro(i).Lbeta(:,:,inds),[],3);
    elseif(strcmp(SL(SLInd).StimHemi,'R'))
        Beta.C(end+1,:,:) = max(SLNeuro(i).Lbeta(:,:,inds),[],3);
        Beta.I(end+1,:,:) = max(SLNeuro(i).Rbeta(:,:,inds),[],3);
    end
   
    if(strcmp(SL(SLInd).Animal,'Igor'))
        StimSites = SL(SLInd).Stim_Loc;
        for s = 1:length(StimSites)
            site = find(cell2mat(cellfun(@(x) strcmp(x,StimSites{s}), SLNeuro(i).chnm,'UniformOutput',0)));
            if(~isempty(site))
                Beta.I(end,:,site) = nan; Beta.C(end,:,site) = nan;
            end
        end
    end    
    StimHemi(end+1) = SL(SLInd).StimHemi;

end

Lchns = 1:(size(Beta.C,3)/2);
Rchns = Lchns + (size(Beta.C,2)/2);
IBetaI = []; CBetaI = [];
for i = 1:size(Beta.C,1)
    if(strcmp(SL(SLInd).StimHemi,'L'))
        IBetaI(end+1,:) = Beta.I(i,2,Lchns)-Beta.I(i,1,Lchns);
        CBetaI(end+1,:) = Beta.I(i,2,Rchns)-Beta.I(i,1,Rchns);
    elseif(strcmp(SL(SLInd).StimHemi,'R'))
        CBetaI(end+1,:) = Beta.I(i,2,Lchns)-Beta.I(i,1,Lchns);
        IBetaI(end+1,:) = Beta.I(i,2,Rchns)-Beta.I(i,1,Rchns);
    end
end

IBetaI = IBetaI(:)'; CBetaI = CBetaI(:)';

[~,j] =hampel(IBetaI);
HIBetaI = IBetaI(~j);
[~,j] =hampel(CBetaI);
HCBetaI = CBetaI(~j);

AllBeta = [HCBetaI,HIBetaI];
groups = [zeros(1,length(HCBetaI)),ones(1,length(HIBetaI))];
figure; boxplot(AllBeta,groups,'labels',{'Contra Trials','Ipsi Trials'}); 
box off; ylabel('\Delta PMBR'); title('Ipsi Hemisphere')
ylim([-5,10])



% RBetaL = Beta.L(:,Rchns); RBetaL = RBetaL(:)';
% RBetaR = Beta.R(:,Rchns); RBetaR = RBetaR(:)';
% LBetaR = Beta.R(:,Lchns); LBetaR = LBetaR(:)';
% LBetaL = Beta.L(:,Lchns); LBetaL = LBetaL(:)';
% 
% [~,j] =hampel(LBetaL);
% HLBetaL = LBetaL(~j);
% [~,j] =hampel(LBetaR);
% HLBetaR = LBetaR(~j);
% [~,j] =hampel(RBetaL);
% HRBetaL = RBetaL(~j);
% [~,j] =hampel(RBetaR);
% HRBetaR = RBetaR(~j);
% 
% AllBeta = [HRBetaR,HRBetaL,HLBetaR,HLBetaL];
% groups = [zeros(1,length(HRBetaR)),ones(1,length(HRBetaL)),2*ones(1,length(HLBetaR)),3*ones(1,length(HLBetaL))];
% figure; boxplot(AllBeta,groups); xticklabels([]); box off; ylabel('\Delta PMBR')








