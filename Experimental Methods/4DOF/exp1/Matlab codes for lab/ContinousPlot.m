function ContinousPlot()
Fs=5e3;T=1;FileName='F1';


%function Dan_Cont_Plot_Data(Fs,T,FileName)
% Continuous data acquiring and plotting
% [Fs] = Hz;        Sampling frequency
% [T] = sec;        Refreshing interval;
% FileName - Must be a string, otherwise will not save

% V2A = 0.4/10;           % Voltage to Amperage constant
% A2N = 7.35;             % Amperage to N
% VCC = V2A*A2N;          % Voltage to N

if Fs > 5e3
    disp('\n*** Input Fs exceeds permitted value ***\n'); 
    return
end

%% Connection
% Create an NI session object and add five analog input channels and one analog output channels
s = daq.createSession('ni');
s.Rate = Fs;                                        % Define sampling frequency
s.DurationInSeconds = Fs*T;
% DevID = daqhwinfo('nidaq','InstalledBoardIds');     % ID Number of NI USB device
% Dev = char(DevID);
d = daq.getDevices;
Dev = d.ID;
s.addAnalogInputChannel(Dev, 0:4, 'Voltage');       % Add five input channels 0-4;->1-5

for kk = 1:4                                        % Optic sensors channels
    s.Channels(kk).InputType='SingleEnded';
	% Min Range is -0.20 to +0.20 Volts
	% Max Range is -10   to 10 Volts
	s.Channels(kk).Range = [-5 5];
end 
s.Channels(5).InputType='SingleEnded';              % Voice Coil Amperage channel
%% Add listener for continuous data acquiring
lh = addlistener(s,'DataAvailable', @ ContPlotData);    % Plot continuous data

fid1 = fopen('data.bin','w');                           % Log Data
lh = s.addlistener('DataAvailable',@(src, event)logData(src, event, fid1));

s.NotifyWhenDataAvailableExceeds = round(Fs*T);
s.IsContinuous = true ;

s.startBackground();

%% Stop Session, plot and save data
while 1
    reply = input('To STOP? Y/N [Y]:','s');
    if strcmpi(reply,'y')==1,
        stop(s);
        delete(lh);
        fclose(fid1);
        fid2 = fopen('data.bin','r');
        [data,~]= fread(fid2,[6,inf],'double');
        fclose(fid2);
        s.stop();
        fclose('all');
        
        % Plot
        time = data(1,:)';
        data = data(2:6,:)';
        delete data.bin
        
        figure(2)
        plot(time,data(:,1:4),'.')
        grid minor
        xlim([0 time(end)])
        title('Sensors reading')
        xlabel('[sec]')
        ylabel('[V]')
        legend('S_1','S_2','S_3','S_4')
        
        figure(3)
        plot(time,data(:,5),'.')
        grid minor
        xlim([0 time(end)])
        title('VC - Force')
        xlabel('[sec]')
        ylabel('[V]')
        
        %% save data
        if ischar(FileName)
            save(FileName, 'data','time');
        end
        break;
    end
end

end

