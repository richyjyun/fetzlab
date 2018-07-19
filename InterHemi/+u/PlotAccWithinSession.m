day = 9;
figure;
window = [-0.4,0.4];
[accs, ts] = a.GetAverageAccRT(UbiSL(day),window);
subplot(1,2,2); % av accel R
line(ts, accs(:,1), 'color', 'k', 'linewidth', 1.5), hold on
line(ts, accs(:,2), 'color', 'b', 'linewidth', 1.5, 'linestyle', ':'), hold on
line(ts, accs(:,3), 'color', 'r', 'linewidth', 1.5, 'linestyle', '--'), hold off
title('Average Ipsi Acceleration')
yticks('')
xlabel('Time (s)')
xlim(window)

subplot(1,2,1) % av accel R
line(ts, accs(:,4), 'color', 'k', 'linewidth', 1.5), hold on
line(ts, accs(:,5), 'color', 'b', 'linewidth', 1.5, 'linestyle', ':'), hold on
line(ts, accs(:,6), 'color', 'r', 'linewidth', 1.5, 'linestyle', '--'), hold off
title('Average Contra Acceleration')
legend({'PRE','STIM','POST'},'Location','Northeast'); legend('boxoff') 
yticks('');
xlabel('Time (s)')
xlim(window)
