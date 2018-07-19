% load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataIgorFinal.mat');
% load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataIgorNeuro.mat');


dates = extractfield(SL,'Date');
Bwindow = [0.3,1]; Gwindow = [-0.5,0.5];
Btime = 0.5; Gtime = 0;
Beta = struct('C',[],'I',[]);
Gamma = struct('C',[],'I',[]);
StimHemi = []; Inds = [];
for i = 1:length(SLNeuro)
    Binds = SLNeuro(i).tneuro > Bwindow(1) & SLNeuro(i).tneuro < Bwindow(2);
    Ginds = SLNeuro(i).tneuro > Gwindow(1) & SLNeuro(i).tneuro < Gwindow(2);
    Bshift = find(Binds,1)-1;
    Gshift = find(Ginds,1)-1;
    
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
    
    for chn = 1:size(SLNeuro(i).Rbeta,2)
        for epoch = 1:size(SLNeuro(i).Rbeta,1)
            if(strcmp(SL(SLInd).StimHemi,'L'))
                
                [~,temp] = findpeaks(squeeze(SLNeuro(i).Rbeta(epoch,chn,:)));
                peaks = SLNeuro(i).tneuro(temp); [~,temp] = min(abs(SLNeuro(i).tneuro(temp)-Btime));
                Beta.C(i,end+1) = peaks(temp);
                
                [~,temp] = findpeaks(squeeze(SLNeuro(i).Lbeta(epoch,chn,:)));
                peaks = SLNeuro(i).tneuro(temp); [~,temp] = min(abs(SLNeuro(i).tneuro(temp)-Btime));
                Beta.I(i,end+1) = peaks(temp);
                
                [~,temp] = findpeaks(squeeze(SLNeuro(i).Rgamma(epoch,chn,:)));
                peaks = SLNeuro(i).tneuro(temp); [~,temp] = min(abs(SLNeuro(i).tneuro(temp)-Gtime));
                Gamma.C(i,end+1) = peaks(temp);
                
                [~,temp] = findpeaks(squeeze(SLNeuro(i).Lgamma(epoch,chn,:)));
                peaks = SLNeuro(i).tneuro(temp); [~,temp] = min(abs(SLNeuro(i).tneuro(temp)-Gtime));
                Gamma.I(i,end+1) = peaks(temp);
                
            elseif(strcmp(SL(SLInd).StimHemi,'R'))
                [~,temp] = findpeaks(squeeze(SLNeuro(i).Lbeta(epoch,chn,:)));
                peaks = SLNeuro(i).tneuro(temp); [~,temp] = min(abs(SLNeuro(i).tneuro(temp)-Btime));
                Beta.C(i,end+1) = peaks(temp);
                
                [~,temp] = findpeaks(squeeze(SLNeuro(i).Rbeta(epoch,chn,:)));
                peaks = SLNeuro(i).tneuro(temp); [~,temp] = min(abs(SLNeuro(i).tneuro(temp)-Btime));
                Beta.I(i,end+1) = peaks(temp);
                
                [~,temp] = findpeaks(squeeze(SLNeuro(i).Lgamma(epoch,chn,:)));
                peaks = SLNeuro(i).tneuro(temp); [~,temp] = min(abs(SLNeuro(i).tneuro(temp)-Gtime));
                Gamma.C(i,end+1) = peaks(temp);
                
                [~,temp] = findpeaks(squeeze(SLNeuro(i).Rgamma(epoch,chn,:)));
                peaks = SLNeuro(i).tneuro(temp); [~,temp] = min(abs(SLNeuro(i).tneuro(temp)-Gtime));
                Gamma.I(i,end+1) = peaks(temp);
            end
            
        end
    end
    
    
    %     if(strcmp(SL(SLInd).StimHemi,'L'))
    %         [~,temp] = max(SLNeuro(i).Rbeta(:,:,Binds),[],3);
    %         Beta.C(end+1,:) = SLNeuro(i).tneuro(temp(:)+Bshift);
    %         [~,temp] = max(SLNeuro(i).Lbeta(:,:,Binds),[],3);
    %         Beta.I(end+1,:) =  SLNeuro(i).tneuro(temp(:)+Bshift);
    %         [~,temp] = max(SLNeuro(i).Rgamma(:,:,Ginds),[],3);
    %         Gamma.C(end+1,:) = SLNeuro(i).tneuro(temp(:)+Gshift);
    %         [~,temp] = max(SLNeuro(i).Lgamma(:,:,Ginds),[],3);
    %         Gamma.I(end+1,:) = SLNeuro(i).tneuro(temp(:)+Gshift);
    %     elseif(strcmp(SL(SLInd).StimHemi,'R'))
    %         [~,temp] = max(SLNeuro(i).Lbeta(:,:,Binds),[],3);
    %         Beta.C(end+1,:) = SLNeuro(i).tneuro(temp(:)+Bshift);
    %         [~,temp] = max(SLNeuro(i).Rbeta(:,:,Binds),[],3);
    %         Beta.I(end+1,:) =  SLNeuro(i).tneuro(temp(:)+Bshift);
    %         [~,temp] = max(SLNeuro(i).Lgamma(:,:,Ginds),[],3);
    %         Gamma.C(end+1,:) = SLNeuro(i).tneuro(temp(:)+Gshift);
    %         [~,temp] = max(SLNeuro(i).Rgamma(:,:,Ginds),[],3);
    %         Gamma.I(end+1,:) = SLNeuro(i).tneuro(temp(:)+Gshift);
    %     end
    
    %     if(strcmp(SL(SLInd).Animal,'Ubi'))
    %         StimSites = strsplit(SL(SLInd).Stim_Loc,'/');
    %         StimSites{2} = [StimSites{1}(1:4),StimSites{2}];
    %     else
    %         StimSites = SL(SLInd).Stim_Loc;
    %     end
    %     for s = 1:length(StimSites)
    %         site = find(cell2mat(cellfun(@(x) strcmp(x,StimSites{s}), SLNeuro(i).chnm,'UniformOutput',0)));
    %         if(~isempty(site))
    %             Beta.I(end,:,site) = nan; Beta.C(end,:,site) = nan;
    %         end
    %     end
    %
    %     StimHemi(end+1) = SL(SLInd).StimHemi;
    %     Inds(end+1) = SLInd;
end

figure;
subplot(2,2,1); histogram(Beta.C(:),80);
subplot(2,2,3); histogram(Beta.I(:),80);
subplot(2,2,2); histogram(Gamma.C(:),80);
subplot(2,2,4); histogram(Gamma.I(:),80);



