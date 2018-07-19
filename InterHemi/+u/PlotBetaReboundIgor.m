% function PlotBetaRebound(SLNeuro,SL,fname)

% fname = fullfile('F:\Dropbox\repos\abogaard\efetz\U\packets', 'BetaGammaUbi.ps');
% if(exist(fname))
%     delete(fname);
% end

% load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataIgorFinal.mat');
% load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataIgorNeuro.mat');

dates = extractfield(SL,'Date'); window = 1.5;
Beta = struct('C',[],'I',[]); 
StimHemi = []; Inds = [];
for i = 1:length(SLNeuro)
    inds = SLNeuro(i).tneuro > 1 & SLNeuro(i).tneuro < window;
    
    %% find proper SL, get change in RT
    % Doing post - pre for igor, can change to cond-pre for ubi
    SLInd = find(strcmp(dates,SLNeuro(i).Date)); trig = SL(SLInd).trig1;
    D = str2num(SL(SLInd).Date);
    if(~strcmp(SL(SLInd).Condition(1:6),'Contra'))
        continue;
    end    
 
    if(size(SLNeuro(i).Rbeta,1)<3 || size(SLNeuro(i).Lbeta,1)<3)
        continue;
    end

    if(strcmp(SL(SLInd).StimHemi,'L'))
        Beta.C(end+1,:,:) = max(SLNeuro(i).Rbeta(:,:,inds),[],3)-min(SLNeuro(i).Rbeta(:,:,inds),[],3);
        Beta.I(end+1,:,:) = max(SLNeuro(i).Lbeta(:,:,inds),[],3)-min(SLNeuro(i).Lbeta(:,:,inds),[],3);
    elseif(strcmp(SL(SLInd).StimHemi,'R'))
        Beta.C(end+1,:,:) = max(SLNeuro(i).Lbeta(:,:,inds),[],3)-min(SLNeuro(i).Lbeta(:,:,inds),[],3);
        Beta.I(end+1,:,:) = max(SLNeuro(i).Rbeta(:,:,inds),[],3)-min(SLNeuro(i).Rbeta(:,:,inds),[],3);
    end
    
    if(strcmp(SL(SLInd).Animal,'Ubi'))
        StimSites = strsplit(SL(SLInd).Stim_Loc,'/');
        StimSites{2} = [StimSites{1}(1:4),StimSites{2}];
    else
        StimSites = SL(SLInd).Stim_Loc;
    end
    for s = 1:length(StimSites)
        site = find(cell2mat(cellfun(@(x) strcmp(x,StimSites{s}), SLNeuro(i).chnm,'UniformOutput',0)));
        if(~isempty(site))
            Beta.I(end,:,site) = nan; Beta.C(end,:,site) = nan;
        end
    end
    
    StimHemi(end+1) = SL(SLInd).StimHemi;
    Inds(end+1) = SLInd;
end

if(strcmp(SL(SLInd).Animal,'Ubi'))
    Rchns = 1:(size(Beta.C,3)/2);
    Lchns = Rchns + Rchns(end);
else
    Lchns = 1:(size(Beta.C,3)/2);
    Rchns = Lchns + Lchns(end);
end
IBetaI = []; CBetaI = []; IBetaC = []; CBetaC = [];
for i = 1:size(Beta.C,1)
    if(strcmp(SL(Inds(i)).StimHemi,'L'))
        IBetaI(end+1,:) = Beta.I(i,2,Lchns)-Beta.I(i,1,Lchns);
        CBetaI(end+1,:) = Beta.C(i,2,Lchns)-Beta.C(i,1,Lchns);
        IBetaC(end+1,:) = Beta.I(i,2,Rchns)-Beta.I(i,1,Rchns);
        CBetaC(end+1,:) = Beta.C(i,2,Rchns)-Beta.C(i,1,Rchns);
    elseif(strcmp(SL(Inds(i)).StimHemi,'R'))
        IBetaI(end+1,:) = Beta.I(i,2,Rchns)-Beta.I(i,1,Rchns);
        CBetaI(end+1,:) = Beta.C(i,2,Rchns)-Beta.C(i,1,Rchns);
        IBetaC(end+1,:) = Beta.I(i,2,Lchns)-Beta.I(i,1,Lchns);
        CBetaC(end+1,:) = Beta.C(i,2,Lchns)-Beta.C(i,1,Lchns);
    end
end

IBetaI = IBetaI(:)'; CBetaI = CBetaI(:)';
IBetaC = IBetaC(:)'; CBetaC = CBetaC(:)';

% [~,j] =hampel(IBetaI);
% HIBetaI = IBetaI(~j);
% [~,j] =hampel(CBetaI);
% HCBetaI = CBetaI(~j);
% [~,j] =hampel(IBetaC);
% HIBetaC = IBetaC(~j);
% [~,j] =hampel(CBetaC);
% HCBetaC = CBetaC(~j);

AllBeta = [CBetaC,IBetaC,CBetaI,IBetaI];
groups = [zeros(1,length(CBetaC)),ones(1,length(IBetaC)),2*ones(1,length(CBetaI)),3*ones(1,length(IBetaI))];
figure; boxplot(AllBeta,groups,'labels',{'ContT ContH','IpsiT ContH','ContT IpsiH','IpsiT IpsiH'}); 
box off; ylabel('\Delta PMBR'); 
% ylim([-5,10])


% figure; plot(SLNeuro(1).tneuro,squeeze(SLNeuro(9).Lbeta(1,15,:)))
% hold on; plot(SLNeuro(1).tneuro,squeeze(SLNeuro(9).Lbeta(2,15,:)))





