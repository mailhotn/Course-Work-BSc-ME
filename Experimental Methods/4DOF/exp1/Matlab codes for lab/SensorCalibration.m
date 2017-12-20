function SensorCalibration()
A=10;F=1;Fs=5e3;N=10;t0=5;
%FileName=['Calibration_SA',regexprep(datestr(now),':','-')];
Answer=listdlg('ListString',{'Mass1','Mass2','Mass3','Mass4'},...
    'PromptString','Which mass?');
FileName=['Calib_',num2str(Answer)];

%function Dan_Sin_Laser_call(F,Fs,A,N,t0,FileName)
% Acquire Data and Generate sin signal simultaneously
% [A] = V;          Amplitude of input excitation;
% [F] = Hz;         Sine frequency;
% [Fs] = Hz;        sampling frequency;
% [N] = #;          Number of cycles
% [t0] = sec        time to wait before the sine signal starts

% FileName - Must be a string, otherwise will not save

FsM = 50e3;               % Max Rate is 50 kHz in the current NI configuration
A = abs(A);

% clear all; close all
AM = 10;     % Max Value - code protection;
if abs(A) > AM
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

s.addAnalogInputChannel(Dev, 0:5, 'Voltage');       % Add five input channels 0-4;->1-5
for kk = 1:4                                        % Optic sensors channels
    s.Channels(kk).InputType='SingleEnded';
	% Min Range is -0.20 to +0.20 Volts
	% Max Range is -10   to 10 Volts
	s.Channels(kk).Range = [-5 5];
end 
s.Channels(5).InputType='SingleEnded';              % Voice Coil Amperage channel
s.Channels(6).InputType='SingleEnded';              % Keyence Laser channel
s.Channels(6).Range = [-5 5];
s.addAnalogOutputChannel(Dev, 0, 'Voltage');        % Voice Coil voltage channel

% Define listener for contionuse plotting
lh = addlistener(s,'DataAvailable', @ContPlotData);

%% Define Output signal
tt = 0:1/Fs:N/F; 
VCSignal  = [zeros(t0*Fs,1); A*sin(2*pi*F*tt.')]; 	% Force excitation
s.queueOutputData(VCSignal);         

s.NotifyWhenDataAvailableExceeds = Fs;                                  % update graph every sec

%% Generate signal and a acquire and store data - Blocks Matlab while doing so
[data, time] = s.startForeground;
delete(lh)

% plot data
figure(2);
plot(time,data(:,[1:4 6]),'.')
grid minor
xlim([0 time(end)])
title(sprintf('Sensors reading, VC amp = %.3f[V] f = %.3f[Hz]',A,F))
xlabel('[sec]')
ylabel('[V]')
legend('M_1','M_2','M_3','M_4','Laser')

figure(3);
plot(time,data(:,5),'.')
grid minor
xlim([0 time(end)])
title('VC - Force')
xlabel('[sec]')
ylabel('[V]')

%% save data
if ischar(FileName)
    save(FileName, 'data','time','A','Fs','F','t0','Answer');
end

%% fit linear curve:
ind = find(time>t0);
t = time(ind);
laser = data(ind,6)-mean(data(ind,6));
Vmass = data(ind,Answer);
figure(4);clf;
subplot 211
plot(time,data(:,6) , t,laser)
legend('Original data','Removed mean','Location','NorthWest')
xlabel('Time [sec]')
ylabel('Displacement [mm]')
grid on

subplot 212
plot(Vmass,laser,'.')
hold on

xlabel('Displacement [mm]')
ylabel('V_mass [V]')
grid on

p = polyfit(Vmass,laser,1);disp(['Fit parameters: ',num2str(p)]);

plot(Vmass,Vmass*p(1)+p(2),'-r');

legend({'V_mass',['fitted curve laser=',num2str(p(1),'%.3f'),'*x+',num2str(p(2),'%.3f')]},'Location','NorthWest')