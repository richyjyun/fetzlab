stLFPs = LFPs.data(:,floor(trialinds(:,1)));


W = sqrtm(inv(cov(stLFPs')));
test1 = W*stLFPs;
figure 
for j = 1:size(LFPs.data,1)
    [c,r,e] = GetWadeChannelPosition(j);
    subplot(10,10,(r-1)*10+c);
    plot(range/fs,test1(j,:),'k'); hold on;
end


figure 
for j = 1:size(LFPs.data,1)
    [c,r,e] = GetWadeChannelPosition(j);
    subplot(10,10,(r-1)*10+c);
    if(j==21)
        plot(range/fs,stLFPs(j,:),'r');
    else
        plot(range/fs,stLFPs(j,:),'k');
    end
%     ylim([-1e-5,1e-5]);
title(num2str(j));
end
%  
% X = bsxfun(@minus,stLFPs',mean(stLFPs'));
% A = X'*X;
% [V,D] = eig(A);
% X = X*V*diag(1./(diag(D)+0.001).^(1/2))*V';

% [E,D] = eig(cov(stLFPs'));
% W = E*inv(sqrtm(D))*E';
% wstLFPs = W*stLFPs;

figure
for j = 1:size(LFPs.data,1)
    [c,r,e] = GetWadeChannelPosition(j);
    subplot(10,10,(r-1)*10+c);
    if(j==21)
        plot(range/fs,wstLFPs(j,:),'r');
    else
        plot(range/fs,wstLFPs(j,:),'k');
    end
    ylim([-0.2,0.2])
end


stLFPs = []; wstLFPs = [];
for j = 1:size(trialinds,2)
    temp = LFPs.data(:,floor(trialinds(:,1)));
    stLFPs(j,:,:) = temp;
    X = temp';
    X = bsxfun(@minus,X,mean(X));
    A = X'*X;
    [V,D] = eig(A);
    X = X*V*diag(1./(diag(D)+0.001).^(1/2))*V';
    
    wstLFPs(j,:,:) = X';
    
end


W = inv(sqrtm(cov(LFPs.data')));
[E,D] = eig(cov(LFPs.data'));
W2 = E*inv(sqrtm(D))*E';

%% Code for whitening plus changing the ratio for plotting

bad = [11,13,15,21,35,66,68,75,77,79,81,95];
% remove bad channels, and average it using those around
temp = LFPs.data;
temp = bpfilt(temp',[15,300],fs,3)';
for j = 1:length(bad)
    [c,r,e] = GetWadeChannelPosition(bad(j));
    avgchns(1) = GetWadeChannel(c-1,r);
    avgchns(2) = GetWadeChannel(c,r-1);
    avgchns(3) = GetWadeChannel(c+1,r);
    avgchns(4) = GetWadeChannel(c,r+1);
    avgchns(avgchns==0) = [];
    [~,ia,~] = intersect(avgchns,bad);
    avgchns(ia) = [];
    
    if(isempty(avgchns))
        temp(bad(j),:) = zeros(1,size(temp,2));
    else
        temp(bad(j),:) = mean(temp(avgchns,:));
    end
    
end

stLFPs = zeros(size(LFPs.data,1),length(range));
% loop through all LFP channels
for j = 1:size(LFPs.data,1)
    
    d = temp(j,:);
    d = d(floor(trialinds));
    d = mean(d,2);
    
    stLFPs(j,:) = d;
    
end

% whiten
[E,D] = eig(cov(temp'));
W = E*inv(sqrtm(D))*E';
wstLFPs = W*stLFPs;

% plot whitened data
figure
for j = 1:96
    [c,r,e] = GetWadeChannelPosition(j);
    subplot(10,10,(r-1)*10+c);
    if(j==21)
        plot(range/fs,wstLFPs(j,:),'r');
    else
        plot(range/fs,wstLFPs(j,:),'k');
    end
    ylim([-0.2,0.2])
end

% find channels with max correlation, find ratio from there
corr = [];
for j = 1:96
    temp = corrcoef(stLFPs(j,:),wstLFPs(j,:));
    corr(j) = abs(temp(1,2));
end
corr(bad) = 0;
ind = find(corr == max(corr));
ratio = abs(max(stLFPs(ind,:))-min(stLFPs(ind,:))) / abs(max(wstLFPs(ind,:))-min(wstLFPs(ind,:)));
wstLFPs = wstLFPs.*ratio;

% plot together
figure;
for j = 1:size(LFPs.data,1)
    [c,r,e] = GetWadeChannelPosition(j);
    subplot(10,10,(r-1)*10+c);
    
    % plot st and wst
    plot(range/fs,stLFPs(j,:),'k'); hold on;
    plot(range/fs,wstLFPs(j,:),'g');
    
    if(j~=trigChns(i))
        yl(j,:) = ylim;
    end
    
end

YLIM = [nanmedian(yl(:,1)),nanmedian(yl(:,2))].*5/4; % get median limits and increase by 20%

% set y lims and other parts
for j = 1:size(LFPs.data,1)
    [c,r,e] = GetWadeChannelPosition(j);
    subplot(10,10,(r-1)*10+c);
    
    ylim(YLIM); hold on;
    line([0 0], YLIM, 'linestyle', ':', 'color', [.7 .7 .7]);
    axis off
    sub_pos = get(gca,'position'); % get subplot axis position
    set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height
    if(j == trigChns(i))
        title(num2str(j),'fontsize',7,'Color','r')
    else
        title(num2str(j),'fontsize',7)
    end
end


