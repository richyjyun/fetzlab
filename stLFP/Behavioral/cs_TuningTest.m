tankpath = 'R:\Fetz Lab\neurowest\spankybackup\OP_DT1_052915\';
blockname = '20150902';

T1 = times(1); T2 = times(2);
Mani = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','Mani'); Mani = Mani.streams.Mani;
Mani.data = [Mani.data(5,:);Mani.data(3,:)]; % 5 is FE, 3 is RU. same as current set up in needing to flip FE.
Snips = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNeu'); Snips = Snips.snips.eNeu;

chn = 58; code = 0;
fs = Mani.fs;

ind = Snips.chan == chn & Snips.sortcode == code;
neu = Snips.data(ind,:);
spk = (Snips.ts(ind)' - T1)*fs;

window = 0.1; dt = 0.02; %time window and steps to look at
[~,avg] = SpatialRate(Mani,spk,fs,window,dt,1);

% avg(avg == 0) = nan;
% for i = 1:size(avg,1)
%     for j = 1:size(avg,2)
%         if(i == 1 || j==1 || i == size(avg,1) || j==size(avg,2))
%             avg(i,j) = nan;
%             continue;
%         else
%             I = [i-1,i+1,i,i];
%             J = [j,j,j+1,j-1];
%             if(sum(isnan(diag(avg(I,J)))) >= 3)
%                 avg(i,j) = nan;
%             end
%         end
%     end
% end

[direction,magnitude] = DirectionalTuning(Mani,spk,fs,window,dt,1);
figure; subplot(2,1,1); histogram(direction);
title('Directional Tuning'); xlabel('Angle (radians)');

edges = -180:6:180; edges = edges*pi/180;

[N,~] = histcounts(direction,edges);
subplot(2,1,2); polarplot(edges,[N,N(1)]);



for i = 1:96
    disp([num2str(i),'-',num2str(sum(Snips.chan==i))]);
end