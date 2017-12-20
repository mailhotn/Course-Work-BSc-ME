function SinInput()
F = 9.55; Fs=5e3; A=5; N=50; DC=0; t0=2; FileName='SinInput_400';
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
figure(1)
plot(time,data(:,1:4),'.')
grid minor
xlim([0 time(end)])
title(sprintf('Sensors reading, VC amp = %.3f[V] dc = %.3f[V] f = %.3f[Hz]',A,DC,F))
xlabel('[sec]')
ylabel('[V]')
legend('M_1','M_2','M_3','M_4')

figure(2)
plot(time,data(:,5),'.')
grid minor
xlim([0 time(end)])
title('VC - Force')
xlabel('[sec]')
ylabel('[V]')

if exist('All_Cs.mat','file')
    load('All_Cs.mat');
    figure(3)
    plot(time,polyval([C3(1) C2(1) C1(1) C0(1)],data(:,1)),...
        time,polyval([C3(2) C2(2) C1(2) C0(2)],data(:,2)),...
        time,polyval([C3(3) C2(3) C1(3) C0(3)],data(:,3)),...
        time,polyval([C3(4) C2(4) C1(4) C0(4)],data(:,4)));
    legend('Mass1','Mass2','Mass3','Mass4');
    
    figure(4);
    Mean=polyval([C3(1) C2(1) C1(1) C0(1)],data(:,1))/4+...
        polyval([C3(2) C2(2) C1(2) C0(2)],data(:,2))/4+...
        polyval([C3(3) C2(3) C1(3) C0(3)],data(:,3))/4+...
        polyval([C3(4) C2(4) C1(4) C0(4)],data(:,4))/4;
    Err=abs(polyval([C3(1) C2(1) C1(1) C0(1)],data(:,1))-Mean)+...
        abs(polyval([C3(2) C2(2) C1(2) C0(2)],data(:,2))-Mean)+...
        abs(polyval([C3(3) C2(3) C1(3) C0(3)],data(:,3))-Mean)+...
        abs(polyval([C3(4) C2(4) C1(4) C0(4)],data(:,4))-Mean);
    plot(time,Err);
    if sum(Err)/length(Err)<.15
        disp('Calibration ok!');
    else
        disp(['Calibration not ok, or ',num2str(F),' is not the first mode']);
    end
end
    
        



%% save data
if ischar(FileName)
    save(FileName, 'data','time','A','Fs','DC','F','Fs');
end

end

