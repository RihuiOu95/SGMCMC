breakpts = [linspace(0,18600,50),linspace(19000,2.2e5,50),linspace(2.2e5,2.5e5,50),linspace(2.5e5,length(logz),30)];
detrend_logz = detrend(logz,2,breakpts);
trend = logz - detrend_logz;
figure(1)
hold on
plot(logz,'b');
%plot(trend,'--r','LineWidth',3);
legend({'Raw Data'});
set(gca,'FontSize',16)
%title("Raw Data and Baseline")
hold off
figure(2)
plot(detrend_logz);
set(gca,'FontSize',16)
%title("Detrended Series")