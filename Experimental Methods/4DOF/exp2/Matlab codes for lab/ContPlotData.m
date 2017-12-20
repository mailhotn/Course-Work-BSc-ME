function ContPlotData(src,event)
    

    if length(event.Data)< 10 
        return
    else
        figure(1)
        subplot(3,1,[1 2])
        plot(event.TimeStamps, event.Data(:,1:4),'.-');
        xlim([min(event.TimeStamps) max(event.TimeStamps)]);
        grid on;
        ylabel('[V]');
        xlabel('[sec]');
        title('Acquired Signals');
        grid on
        legend('s1','s2','s3','s4');
        shg

        subplot(3,1,3)    
        plot(event.TimeStamps, event.Data(:,5),'.-');
        xlim([min(event.TimeStamps) max(event.TimeStamps)]);
        grid on;
        ylabel('[V]');
        xlabel('[sec]');
        title('Acquired Signal');
        grid on
        legend('VC current - Force');
    end
end