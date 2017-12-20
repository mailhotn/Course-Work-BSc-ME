function Stepped_Sine_Bode_Nyquist()
Fs=5e3;A=0.5;F1=5;F2=25;df=1.0;N=100;NC=50;Nhar=3;
FileName='bode_FIRST_ZERO(13Hz)';
C0=[0 0 0 0]';C1=[1 1 1 1]';C2=[0 0 0 0];C3=[0 0 0 0 ];
if exist('All_Cs.mat','file')
    load('All_Cs.mat');
    disp('Calibration parameters loaded successfully.');
end
%% Bode Plot using Step Sine

% [Fs] = Hz;        Sampling frequency
% [A] = V;          Amplitude of input excitation Max Valu = 0.5
% [F1] = Hz;        Minmum Frequency frequency
% [F2] = Hz;        Maximum Frequency frequency
% [df] = Hz;        Frequency step
% [N] = #;          Number of cycles at each frequency
% [NC] = #;         Number of cycles which are assumed as transient
% [C0,C1,C2,C3] = mm/V Are the callibration constants of the model: X =
%                   C0+C1*V+C2*V^2+C3*V^3, they should be insetet as vectors
%                   with 4 cells. for example: C1 = [10 10.5 8 3].';
% [Nhar] = #;       Nuber of harmonies to fit
%
% FileName - Must be a string, otherwise will not save

% Fs = 5e3;       % Max Rate is 50 kHz in the current NI configuration

AM = 8;     % Max Value - code protection;
if abs(A)> AM
    disp('\n*** Inpute Voltage exceeds permitted value ***\n'); 
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
s.Channels(5).InputType='SingleEnded';              % Voice Coile Amperage channel
s.addAnalogOutputChannel(Dev, 0, 'Voltage');        % Voice Coile voltage channel

FF = F1:df:F2;                                      % Define frequencies vector

AmpSR = zeros(4,length(FF));
PhaseS = AmpSR;
H = AmpSR;
for ll = 1:length(FF)
    F = FF(ll);
    [AmpSR(:,ll),PhaseS(:,ll),H(:,ll)] = SineGen_LS(s,A,F,Fs,N,NC,Nhar,C0,C1,C2,C3);

    figure(6)
    subplot(2,1,1)
    plot(FF(1:ll),AmpSR(:,1:ll)','--o')
    xlabel('$f$[Hz]','interpreter','latex','FontSize',14)
    ylabel('Amplitude Ratio [mm/N]','interpreter','latex','FontSize',14)
    title('Bode','interpreter','latex','FontSize',14)
    xlim([FF(1) FF(end)]);
    grid on;
    hl = legend('$m_1$','$m_2$','$m_3$','$m_4$');
    set(hl,'interpreter','latex','FontSize',14);
    subplot(2,1,2)
    plot(FF(1:ll),unwrap(PhaseS(:,1:ll)')/pi*180,'--o')
    xlabel('$f$[Hz]','interpreter','latex','FontSize',14)
    ylabel('Phase [deg]','interpreter','latex','FontSize',14)
    xlim([FF(1) FF(end)]);
    grid on;
    hl = legend('$m_1$','$m_2$','$m_3$','$m_4$');
    set(hl,'interpreter','latex','FontSize',14);
end

AmpSR = AmpSR.';
PhaseS = PhaseS.';
H = H.';

%% plot the frequency response as a Nyquist plot
for kk = 1:4
    figure
    plot(H(:,kk),'--o','LineWidth',2)
    xlabel('Real${H}$','interpreter','latex','FontSize',14)
    ylabel('Imaginary${H}$','interpreter','latex','FontSize',14)
    title(sprintf('M No. %d',kk),'interpreter','latex','FontSize',14)
    axis equal
end

%% save data
if ischar(FileName)
    save(FileName);
end

end

