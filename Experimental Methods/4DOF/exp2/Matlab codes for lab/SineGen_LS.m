function [AmpR,Phase,H] = SineGen_LS(s,A,F,Fs,N,NC,Nhar,C0,C1,C2,C3)
% Generate Sine Signal and Acquire data from the VC and the optic sensors.
% Computing amplitudes and phases using linear LS

V2A = 0.4/10;           % Voltage to Amperage constant
A2N = 7.35;             % Amperage to N
VCC = V2A*A2N;          % Voltage to N

AmpR = zeros(4,1);
Phase = AmpR;
H = AmpR;
Amps = zeros(Nhar,1);
%% Define Output signal
tt = 0:1/Fs:N/F; 
VCSignal  = A*sin(2*pi*F*tt.');                % Force excitaion
s.queueOutputData(VCSignal);         

%% Generate signal and a acquire and store data - Blocks Matlab while duing so
[data, time] = s.startForeground;
% data = dtrend(data);                % remove trend
data(time<NC/F,:) = [];
tt(time<NC/F) = [];

%% Build the model
B = [cos(F*(2*pi)*tt.'), sin(F*(2*pi)*tt.')];
for jj = 2:Nhar
    B = [B, cos(jj*F*(2*pi)*tt.'), sin(jj*F*(2*pi)*tt.')];
end
data(:,5) = dtrend(data(:,5)) * VCC;
OutPut = B\data(:,5);
figure(1)
subplot(1,3,[1 2])
h = plot(tt, [data(:,5), B*OutPut, data(:,5)-B*OutPut]);
set(h(1),'Color','k');
set(h(2),'Color','b')
set(h(3),'Color','r');
legend(h,'Messurments','Fitted','Error');
title('Input Force')
xlabel('t[sec]')
ylabel('F[N]');
xlim([tt(1) tt(end)]);
grid on
subplot(1,3,3)

for jj = 1:Nhar
    Amps(jj) = sqrt(OutPut(2*jj-1)^2+OutPut(2*jj)^2);
end
bar(1:Nhar,Amps)

for kk = 1:4
    MesuredSignal = C0(kk)+C1(kk)*data(:,kk)+C2(kk)*data(:,kk).^2+C3(kk)*data(:,kk).^3;
    a  = B\MesuredSignal;
    figure(kk+1)
    subplot(1,3,[1 2])
    h = plot(tt, [MesuredSignal, B*a, MesuredSignal-B*a]);
    set(h(1),'Color','k');
    set(h(2),'Color','b')
    set(h(3),'Color','r');
    legend(h,'Messurments','Fitted','Error');
    title(sprintf('M No. %d',kk))
    xlabel('t[sec]')
    ylabel('X[mm]');
    grid on
    xlim([tt(1) tt(end)]);
    for jj = 1:Nhar
        Amps(jj) = sqrt(a(2*jj-1)^2+a(2*jj)^2);
    end
    subplot(1,3,3)
    bar(1:Nhar,Amps)
    shg

    %% Extract data
    AmpR(kk) = sqrt(a(1)^2+a(2)^2)/(A*VCC);
    Phase(kk) = angle((a(1)-1i*a(2))/(OutPut(1)-1i*OutPut(2)));
    H(kk) = (a(1)-1i*a(2))/(OutPut(1)-1i*OutPut(2));
end

end

