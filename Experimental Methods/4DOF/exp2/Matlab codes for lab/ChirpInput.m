function ChirpInput()
Fs=5e3;A=5;f0=3;f1=50;T=10;FileName='Ch3_50_60sec';


%function Dan_Chirp(Fs,A,f0,f1,T,FileName)
% Acquire Data and Generate chirp signal simultaneously

% [Fs] = Hz;        Sampling frequency
% [A] = V;       	Amplitude of input excitation Max Value = 10;
% [f0] = Hz;      	Chirp, minimum ftequency
% [f1] = Hz         Chirp, maximum ftequency
% [T] = sec;        Duration of the signal
% FileName - Must be a string, otherwise will not save

% Fs = 5e3;               % Max Rate is 50 kHz in the current NI configuration
% V2A = 0.4/10;           % Voltage to Amperage constant
% A2N = 7.35;             % Amperage to N
% VCC = V2A*A2N;          % Voltage to N

% clear all; close all

AM = 10;     % Max Value - code protection;
if abs(A) > AM
    disp('\n*** Input Voltage exceeds permitted value ***\n'); 
    return
end

Fs = abs(Fs);
if Fs> 50e3
    disp('\n*** Inpute sampling frequency exceeds permitted value ***\n'); 
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

%% Define listener for continuous plotting
lh = addlistener(s,'DataAvailable', @ContPlotData);

%% Define Output signal
N = T*Fs;
u = chirp(linspace(0,T,N)',f0,T,f1);

%VCSignal  = [zeros(2*Fs,1) ; A*ones(T1*Fs,1); zeros(T2*Fs,1)];	% Force excitation
s.queueOutputData(A*u);
s.NotifyWhenDataAvailableExceeds = round(Fs/2);           % update graph twice a second

%% Generate signal and a acquire and store data - Blocks Matlab while doing so
[data, time] = s.startForeground;
delete(lh)

% plot data
figure
plot(time,data(:,1:4),'.')
grid minor
xlim([0 time(end)])
title(sprintf('Sensors reading, VC amp = %.3f[V]',A))
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
    save(FileName, 'data','time','Fs','A','f0','f1','T','u');
end

end

