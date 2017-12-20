function SinInput()
% F =10.1; A=5;
Fs=5e3;  N=50; DC=0; t0=4; 
Answer=inputdlg({'Please enter actuation frequency:','Amplitude'});
F=str2double(Answer{1});
A=str2double(Answer{2});
FileName=regexprep(['SinInput_',num2str(F)],'\.','_');
if abs(F-10)<2 && A>3
    error('Ampiltude too high for the first eigen frequency');
end
    
% Acquire Data and Generate sin signal simultaneously
% [A] = V;          Amplitude of input excitation;
% [DC] = V;         Additional DC to the signal;
% [F] = Hz;         Sine frequency;
% [Fs] = Hz;        sampling frequency;
% [N] = #;          Number of cycles
% [t0] = sec        time to wait after signal ends

% FileName - Must be a string, otherwise will not save

FsM = 50e3;               % Max Rate is 50 kHz in the current NI configuration
A = abs(A);
% V2A = 0.4/10;           % Voltage to Amperage constant
% A2N = 7.35;             % Amperage to N
% VCC = V2A*A2N;          % Voltage to N

% clear all; close all
AM = 10;     % Max Value - code protection;
if abs(A)+abs(DC) > AM
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
tt = 0:1/Fs:N/F; 
VCSignal  = [zeros(2*Fs,1); DC + A*sin(2*pi*F*tt.'); zeros(t0*Fs,1)]; 	% Force excitation
s.queueOutputData(VCSignal);         

s.NotifyWhenDataAvailableExceeds = Fs;                                  % update graph every sec

%% Generate signal and a acquire and store data - Blocks Matlab while doing so
[data, time] = s.startForeground;
delete(lh)

% plot data
figure
plot(time,data(:,1:4),'.')
grid minor
xlim([0 time(end)])
title(sprintf('Sensors reading, VC amp = %.3f[V] dc = %.3f[V] f = %.3f[Hz]',A,DC,F))
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
    save(FileName, 'data','time','A','Fs','DC','F','Fs');
end

end

