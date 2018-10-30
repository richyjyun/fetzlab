clear;

%% Load in proper times
tankpath = 'Y:\~NeuroWest\Spanky\SpikeTrigger-180122-105223\';

Change = 'Spanky-180427-134343';


% Change = ['Spanky-180122-130528';
%     'Spanky-180129-131027';
%     'Spanky-180204-155626';
%     'Spanky-180208-133323';
%     'Spanky-180214-140216';
%     'Spanky-180302-162140';
%     'Spanky-180304-154857';
%     'Spanky-180309-132433'];

Pre = cell(1,size(Change,1)); Post = Pre;
for i = 1:8
    blockname = Change(i,:);
    TT = TDT2mat([tankpath,blockname],'TYPE',2);
    Dscm = TT.epocs.Dscm;
    [val,ind] = findpeaks(Dscm.data);
    ind = ind(val>1000); val = val(val>1000);
    tests = [1,3];
    times = [ind(tests)-val(tests),ind(tests)-val(tests)+val(1)];  % times = [ind(1) - val(1),ind(1)-val(1)+15000;ind(3)-val(3),ind(3)-val(3)+10000];
    times = Dscm.onset(times);
    
    [Pre{i},Post{i}] = getCycleTrigLFP(tankpath,blockname,79,times);
end

save('F:\S\Code\u\CycleTrigLFP.mat','Pre','Post','-v7.3');

%%
figure; colors = get(gca,'ColorOrder');
chns = [64]; phs = [];
str1 = ''; str2 = ''; Leg = [];
preStim = []; postStim = [];
idx = [];
for chn = 1:length(chns)
    subplot(2,1,1);
    t = linspace(-100,100,size(Pre{1}{1},1));
    data = cellfun(@(x) x(chns(chn)),Pre);
    data = mean(cat(2,data{:}),2);
    data = data*1e6;
    phase = angle(hilbert(data));
    if(chn==1)
        preStim = phase;
    end
    trough = find(diff(phase) <0);
    [~, index] = min(abs(trough-length(phase)/2));
    %     str1 = sprintf('%s%d Phase @0: %f\n',str1,chns(chn),phase(ceil(length(phase)/2)));
    phs(chn,1) = phase(ceil(length(phase)/2));
    Leg(chn) = plot(t,data,'linewidth',2,'color',colors(chn,:)); hold on;
    yl = ylim;%[-40,40];
    plot([t(trough(index)),t(trough(index))],yl,'--','color',colors(chn,:));
    idx(chn,1) = trough(index);

%     subplot(2,1,2);
%     data = cellfun(@(x) x(chns(chn)),Post);
%     data = mean(cat(2,data{:}),2);
%     data = data*1e6;
%     phase = angle(hilbert(data));
%     if(chn==1)
%         postStim = phase;
%     end
%     trough = find(diff(phase) <0);
%     [~, index] = min(abs(trough-length(phase)/2));
%     %     str2 = sprintf('%s%d Phase @0: %f\n',str2,chns(chn),phase(ceil(length(phase)/2)));
%     phs(chn,2) = phase(ceil(length(phase)/2));
%     plot(t,data,'linewidth',2,'color',colors(chn,:)); hold on;
%     plot([t(trough(index)),t(trough(index))],yl,'--','color',colors(chn,:));
%     idx(chn,2) = trough(index);

end
% 
% dPhase(1,1) = preStim(idx(1,1)) - preStim(idx(2,1));
% dPhase(1,2) = postStim(idx(1,2)) - postStim(idx(2,2));
% dPhase(2,1) = preStim(idx(1,1)) - preStim(idx(3,1));
% dPhase(2,2) = postStim(idx(1,2)) - (postStim(idx(3,2))+2*pi);
% 
% str1 = sprintf('Trig %sPhase: %0.2f\nCont %sPhase: %0.2f','\Delta',dPhase(1,1),'\Delta',dPhase(2,1));
% 
% subplot(2,1,1); title('Pre');
% legend(Leg,'Stim','Trig','Cont'); legend boxoff;
% ylim(yl); xlim([-15,35]); xlabel('Time (ms)');
% ylabel('Beta (\muV)');
% xl = xlim; yl = ylim;
% text(xl(2),yl(1),str1,'horizontalalignment','right','verticalalignment','bottom');
% 
% str2 = sprintf('Trig %sPhase: %0.2f\nCont %sPhase: %0.2f','\Delta',dPhase(1,2),'\Delta',dPhase(2,2));
% 
% subplot(2,1,2); title('Post');
% ylim(yl); xlim([-15,35]); xlabel('Time (ms)');
% ylabel('Beta (\muV)');
% xl = xlim; yl = ylim;
% text(xl(2),yl(1),str2,'horizontalalignment','right','verticalalignment','bottom');

% subplot(2,2,3);
% xlim([-15,35]); xlabel('Time (ms)'); ylabel('Phase (radians)');
% subplot(2,2,4);
% xlim([-15,35]); xlabel('Time (ms)'); ylabel('Phase (radians)');


% figure;
% subplot(2,1,1);
% t = linspace(-100,100,size(PreStim{1},1));
% data = mean(cat(2,PreStim{:}),2);
% phase{1,1} = angle(hilbert(data));
% str = sprintf('Stim Phase @0: %f\n',phase{1,1}(ceil(length(data)/2)));
% plot(t,data,'linewidth',2); hold on;
%
% data = mean(cat(2,PreTrig{:}),2);
% phase{1,1} = angle(hilbert(data));
% str = sprintf('%sTrig Phase @0: %f\n',str,phase{1,2}(ceil(length(data)/2)));
% plot(t,data,'linewidth',2);
% title('Pre'); legend('Stim','Trig');
% xlabel('Time (ms)'); ylabel('Beta (\muV)');
% ylim([-3e-5,3e-5])
%
% xl = xlim; yl = ylim;
% text(xl(1)+diff(xl)*0.7,yl(1)+diff(yl)*0.1,str);
%
%
% subplot(2,1,2);
% data = mean(cat(2,PreStim{:}),2);
% phase{2,1} = angle(hilbert(data));
% str = sprintf('Stim Phase @0: %f\n',phase{2,1}(ceil(length(data)/2)));
% plot(t,data,'linewidth',2); hold on;
%
% data = mean(cat(2,PostTrig{:}),2);
% phase{2,2} = angle(hilbert(data));
% str = sprintf('%sTrig Phase @0: %f',str,phase{2,2}(ceil(length(data)/2)));
% plot(t,data,'linewidth',2);
% title('Post');
% xlabel('Time (ms)'); ylabel('Beta (\muV)');
% ylim([-3e-5,3e-5])
%
% xl = xlim; yl = ylim;
% text(xl(1)+diff(xl)*0.7,yl(1)+diff(yl)*0.1,str);



