%% Get Calibration
calib = zeros(6,4);
for ii = 1:4
    load(['Calib_' num2str(ii) '.mat']);
    ind = find(time>t0);
    t = time(ind);
    laser = data(ind,6)-mean(data(ind,6));
    Vmass = data(ind,ii);
    [calib(ii,:),S] = polyfit(Vmass,laser,3);
end
calib(5,:) = [0 0 4e-3*7.35 0];
% plot example calibration curve
% load('Calib_3.mat');
% ind = find(time>t0);
% t = time(ind);
% laser = data(ind,6)-mean(data(ind,6));
% Vmass = data(ind,3);
% figure(1)
% plot(Vmass,laser,'.',Vmass,polyval(calib(3,:),Vmass),'Linewidth',2)
% ylabel('Displacement [mm]')
% xlabel('V_ [V]')
% grid on
% legend('V_{mass}','fitted curve','Location','Best')
%% Q2 Sampling Rate
load('SinDownSamp_15_Sampling1000.mat');  % <----- Change
dist = zeros(size(data));
for ii = 1:5
    dist(:,ii) = polyval(calib(ii,:),data(:,ii));
end
data = dist;
figure(1);
plot(time,data(:,1:4),'--o')
grid minor
xlim([0 time(end)])
% title('Down Sampling: F_s = 11 Hz')
xlabel('Time [sec]')
ylabel('Horizontal Displacement[mm]')
legend('M_1','M_2','M_3','M_4')
%% Q2 Fourier
Fs = 60; % <----- Change
time(time<=6) = 0;
time(time>=14) = 0;
ind = find(time);
[s,f] = FourierSeriesDFT(data(ind,1),Fs);
figure(2)
plot(f,abs(s))
xlabel('f [Hz]')
ylabel('Amplitude [mm]')

%% Q2 Regression
Fsample = [1000 60 27 11];
AliasFreq = [15 15 12 4];
for ii = 1:4
    load(['SinDownSamp_15_Sampling' num2str(Fsample(ii)) '.mat']);
    dist = zeros(size(data));
    for jj = 1:5
        dist(:,jj) = polyval(calib(jj,:),data(:,jj));
    end
    data = dist;
    time(time<=6) = 0;
    time(time>=14) = 0;
    ind = find(time);
    data(ind,1) = data(ind,1) - mean(data(ind,1));
    fitObj = fitnlm(time(ind), data(ind,1), @(a,x) ...
        a(1)*sin(2*pi*a(2)*x + a(3)), [0.3 AliasFreq(ii), pi]);
    figure()
    t = linspace(6,14,1000*(14 - 6)); % 1000 samples/sec
    h = plot(t, predict(fitObj,t.'), time(ind), data(ind,1), '--o');
    legend('Least Square Prediction','Raw Data');
    xlabel('Time [Sec]'); ylabel('Amplitude [mm]');
end

%% Q5 Linearity Test 
%%% NOTE: Please run Get Calibration section for this section to work %%%
load('CombOut');
y_comb = data(:,1) - mean(data(:,1));
load('y1Out');
y1 = data(:,1) - mean(data(:,1)); 
load('y2Out');
y2 = data(:,1) - mean(data(:,1));
a1 = 5; a2 = 4; f1 = 5; f2 = 3;
err = y_comb - (a1*y1 + a2*y2);  
PositionsMat = [y1, y2, y_comb, err];
for ii = 1:4
    CalibPositionsMat(:,ii) = polyval(calib(1,:), PositionsMat(:,ii)); %#ok
end
% clean out the bias
Caliby1     = CalibPositionsMat(:,1) - mean(CalibPositionsMat(:,1));
Caliby2     = CalibPositionsMat(:,2) - mean(CalibPositionsMat(:,2));
Caliby_comb = CalibPositionsMat(:,3) - mean(CalibPositionsMat(:,3));
Caliberr    = CalibPositionsMat(:,4) - mean(CalibPositionsMat(:,4));

figure(1);
subplot(2,1,1)
plot(time, Caliby_comb, time, a1*Caliby1 + a2*Caliby2)
xlabel('Time [sec]');
ylabel('Position of mass #1 [mm]');
legend('y_{comb}','a1*y1+a2*y2');
subplot(2,1,2)
plot(time,Caliberr,'r')
xlabel('Time [sec]');
ylabel('Error');
%% Q7 Sine Input
Fs=5e3;
Fr = [9.937 19.62 33.12 42.3];
phasemut = [];
phaseio = [];
for jj = 1:4
    dec = num2str(Fr(jj) - floor(Fr(jj)));
    load(['SinInput_' num2str(floor(Fr(jj))) '_' dec(3:end) '.mat'])
    dist = zeros(size(data));
    for ii = 1:5
        dist(:,ii) = polyval(calib(ii,:),data(:,ii));
    end
    data = dist;
%     % plot data
%     figure(1)
%     plot(time,data(:,1:4),'.')
%     grid minor
%     xlim([0 time(end)])
%     xlabel('[sec]')
%     ylabel('[V]')
%     legend('M_1','M_2','M_3','M_4')
%     
%     figure(2)
%     plot(time,data(:,5),'.')
%     grid minor
%     xlim([0 time(end)])
%     title('VC - Force')
%     xlabel('[sec]')
%     ylabel('Input Force [N]')
    
    time(time<=8) = 0;
    time(time>=11) = 0;
    ind = find(time);
    datas = [];
    for ii = 1:5
        [s,f] = FourierSeriesDFT(data(ind,ii),Fs);
        datas(:,ii) = s.';
    end
    [~, frind] = min(abs(f-Fr(jj)));
    % figure(2)
    % plot(f,abs(s))
    phaseio(jj,:) = angle(datas(frind,1:4)./datas(frind,5));
    phasemut(jj,:) = angle(datas(frind,1:4)./datas(frind,1));
end
phaseio*180/pi
phasemut*180/pi
%% Q9 Model Identification
Fs=5e3;
Fr = 9.937;
load('SinInput_9_937.mat')
    dist = zeros(size(data));
    for ii = 1:5
        dist(:,ii) = polyval(calib(ii,:),data(:,ii));
    end
    data = dist;
mean(damping_coeff(time,data(:,1:4),Fr))
close all
Fr = 19.62;
load('SinInput_19_62.mat')
    dist = zeros(size(data));
    for ii = 1:5
        dist(:,ii) = polyval(calib(ii,:),data(:,ii));
    end
    data = dist;
mean(damping_coeff(time,data(:,1:4),Fr))
close all
Fr = 33.12;
load('SinInput_33_12.mat')
    dist = zeros(size(data));
    for ii = 1:5
        dist(:,ii) = polyval(calib(ii,:),data(:,ii));
    end
    data = dist;
mean(damping_coeff(time,data(:,1:4),Fr))
close all
Fr = 42.3;
load('SinInput_42_3.mat')
    dist = zeros(size(data));
    for ii = 1:5
        dist(:,ii) = polyval(calib(ii,:),data(:,ii));
    end
    data = dist;
mean(damping_coeff(time,data(:,1:4),Fr))
close all