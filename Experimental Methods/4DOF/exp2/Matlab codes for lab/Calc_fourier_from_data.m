function Calc_fourier_from_data(data,time)
% Calculates DFT to input data.
% [data] is a matrix which 4 first columns are the sensors readings.
% [time] is the time vector of the sampling

% calculate sampling frequency
dt = time(2)-time(1); Fs=1/dt;

%% choose data range  to be fitted
figure
plot(time,data(:,1:4))
title('choose 2 points on a single exp. decay part','Interpreter','Latex','FontSize',18)
xlabel('time[sec]','Interpreter','Latex','FontSize',18)
set(gca,'FontSize',14)
set(gcf,'Color','w')
shg
[xx,~] = ginput(2);

xx = sort(xx);
i1 = find(time>=xx(1) & time<=xx(2));

plot(time(i1),data(i1,1:4))
shg
t = time(i1);

for q = 1:4
    % Load data from from a sensor
    y = detrend(data(:,q));
    yt = y(i1);
    
    [Y,f] = FourierSeriesDFT(yt.*hanning(length(yt)),Fs);

    Ydb = 20*log10(abs(Y));

    pks = find(Ydb>-66);
    
    figure
    h = plot(f,Ydb,'--o',f(pks),Ydb(pks),'x');
    set(h(1),'Color',lines(1),'MarkerFaceColor',lines(1),'MarkerEdgeColor',lines(1))
    set(h(2),'MarkerFaceColor','r','MarkerSize',10,'LineWidth',3)
    xlabel('$f$[Hz]','Interpreter','Latex','FontSize',18)
    ylabel('$|Y(\omega)|$[dB]','Interpreter','Latex','FontSize',18)
    title(sprintf('Measurements of mass %d',q),'Interpreter','Latex','FontSize',18)
    if Fs<200
        xlim([1 Fs/2])
    else
        xlim([1 60])
    end
    set(gca,'FontSize',14)
    set(gcf,'Color','w')
    shg
end
end

