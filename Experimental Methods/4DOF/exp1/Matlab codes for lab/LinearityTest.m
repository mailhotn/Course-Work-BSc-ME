function LinearityTest(Fs,VCSignal,FileName)
% Acquire Data and Generate sin signal simultaneously
% [Fs] = Hz;        sampling frequency;
% [N] = #;          Number of cycles
% [VCSignal] = V    Input signal

% FileName - Must be a string, otherwise will not save

FsM = 50e3;               % Max Rate is 50 kHz in the current NI configuration
% V2A = 0.4/10;           % Voltage to Amperage constant
% A2N = 7.35;             % Amperage to N
% VCC = V2A*A2N;          % Voltage to N

% clear all; close all
AM = 10;     % Max Value - code protection;
if max(abs(VCSignal)) > AM
    disp('*** Input Voltage exceeds permitted value ***'); 
    return
end

if abs(Fs) > FsM || Fs<=0
    disp('*** Input Fs value - illegal ***'); 
    return
end

%% Connection
% Create an NI session object and add five analog input channels and one analog output channels
s = daq.createSession('ni');
s.Rate = Fs;                                        % Define sampling frequency
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
s.addAnalogOutputChannel(Dev, 0, 'Voltage');        % Voice Coil voltage channel

% Define listener for contionuse plotting
lh = addlistener(s,'DataAvailable', @ContPlotData);

%% Define Output signal
s.queueOutputData(VCSignal);         
s.NotifyWhenDataAvailableExceeds = Fs;  % update graph every sec

%% Generate signal and a acquire and store data - Blocks Matlab while doing so
[data, time] = s.startForeground;
delete(lh)

% plot data
figure
plot(time,data(:,1:4),'.')
grid minor
xlim([0 time(end)])
xlabel('[sec]')
ylabel('[V]')
legend('M_1','M_2','M_3','M_4')

figure
plot(time,data(:,5),'.')
grid minor
xlim([0 time(end)])
title('VC - Force')
xlabel('[sec]')
ylabel('[V]')

%% save data
if ischar(FileName)
    save(FileName, 'data','time','VCSignal','Fs');
end

end

