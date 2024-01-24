data = load('data/sleep_c3a2.mat'); 
y = data.C3A2;

%%
figure(1);

subplot(1,2,1);
plot(y);
xlabel('Observations (200 Hz)');
ylabel('EEG signal (nV)');
title('Raw Sleep Cycle Data');
xlim([1,4 * 10^6]);
ax = gca;
ax.FontSize = 26;

subplot(1,2,2);
trans_y = log10(y + 1000.0);
plot(trans_y)
xlabel('Observations (200 Hz)');
ylabel('Log EEG signal');
title('Transformed Sleep Cycle Data');
xlim([1,4 * 10^6]);
ax = gca;
ax.FontSize = 26;