function SquareInput()
T=10;A=10;N=4;DC=0;FileName='ExpSq1';

% Acquire Data and Generate square signal simultaneously
% [A] = V;          Amplitude of input excitation;
% [DC] = V;         Additional DC to the signal;
% [T] = sec;        1 Cycle;
% [N] = #;          Number of cycles
% FileName - Must be a string, otherwise will not save

Fs = 5e3;       % Max Rate is 50 kHz in the current NI configuration
% V2A = 0.4/10;           % Voltage to Amperage constant
% A2N = 7.35;             % Amperage to N
% VCC = V2A*A2N;          % Voltage to N

% clear all; close al`l
AM = 10;     % Max Value - code protection;
if abs(A)+abs(DC) > AM
    disp('*** Input Voltage exceeds permitted value ***'); 
    return
end

%% Connection
% Create an NI session object and add five analogue input channels and one analogue output channels
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
tt = 0:1/Fs:N*T;                                        % time vector
VCSignal  = DC + 0.5*A*(1 + square(2*pi/T*tt.'));      	% Force excitation
s.queueOutputData(VCSignal);        

s.NotifyWhenDataAvailableExceeds = T/2*Fs;              % update graph every half cycle

%% Generate signal and a acquire and store data - Blocks Matlab while doing so
[data, time] = s.startForeground;
delete(lh)

% plot data
figure
plot(time,data(:,1:4),'.')
grid minor
xlim([0 time(end)])
title(sprintf('Sensors reading, VC amp = %.3f[V] DC = %.3f[V]',A,DC))
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
    save(FileName, 'data','time','A','T','DC');
end

end

