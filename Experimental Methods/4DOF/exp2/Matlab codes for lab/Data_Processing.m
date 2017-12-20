%% Data Calibration
calib = zeros(5,4);
for ii = 1:4
    load(['Cal12_' num2str(ii) '.mat']);
    ind = find(time>t0);
    t = time(ind);
    laser = data(ind,6) - mean(data(ind,6));
    Vmass = data(ind,ii);
    [calib(ii,:),S] = polyfit(Vmass,laser,3);
    figure()
    subplot(2,1,1)
    plot(t, laser, t, Vmass - mean(Vmass))
    legend('Calibration Sensor [mm]', 'PSD Sensor [V]');
    subplot(2,1,2)
    plot(Vmass,laser, '.', Vmass, polyval(calib(ii,:), Vmass), 'Linewidth', 2)
    ylabel('Displacement [mm]')
    xlabel('V_ [V]')
    grid on
    legend('V_{mass}','fitted curve','Location','Best')
end
calib(5,:) = [0 0 4e-3*7.35 0];

%% Distance-Damping coeff Model
% NOTE: Run Data Calibration section first
close all
Dist = [0 50 100 150 200 250 300].';
for ii = 1:9
    load(['SinInput_' num2str(Dist(ii)) '.mat']);
    for jj = 1:5
        OutSig(:,jj) = polyval(calib(jj,:),data(:,jj)); %#ok
    end
    Zeta(ii,:) = damping_coeff(time, data, F).'; %#ok
end
close all
% load('Zeta.mat');
[FitObj, Stats, Output] = fit(Dist,Zeta,'poly1');
plot(Dist,Zeta,'ro',Dist,feval(FitObj,Dist));
xlabel('Distance [\mum]'); ylabel('\zeta');
fitlineSTR = sprintf(['Linear Fit\n y=' num2str(FitObj.p1)...
    'x ' num2str(FitObj.p2) '\nR^{2}=' num2str(Stats.rsquare)]);
legend('Measured Data', fitlineSTR);
hold on 
errorbarxy(Dist, Zeta, 10, [],{'o','r','r'});
hold off

%% Chirp Response
% NOTE: Run Data Calibration section first
close all
ChirpTime = [5 20 40 60].';
for ii = 1:4
    load(['Ch3_50_' num2str(ChirpTime(ii)) 'sec.mat']);
    for jj = 1:5
        OutSig(:,jj) = polyval(calib(jj,:),data(:,jj));
    end
    figure();
    plot(time, OutSig(:,1:4));
    xlabel('Time [sec]'); ylabel('Displacement [mm]');
    legend('M1','M2','M3','M4');
    clear('OutSig');
end

%% Impulse Input
% NOTE: Run Data Calibration section first
close all
load('ImpulseData.mat');
for jj = 1:5
    OutSig(:,jj) = polyval(calib(jj,:),data(:,jj));
end
plot(time,OutSig);
xlabel('Time [sec]'); ylabel('Displacement [mm]');
legend('M1','M2','M3','M4');
Calc_fourier_from_data(OutSig,time);


%% Q3 Bode, Q4 Nyquist
load ('bode_FIRST_ZERO(13Hz)5_to_25_rough.mat');
f = FF;
H1 = H;
load ('bode_FIRST_ZERO(13Hz)25_to_45_rough.mat');
f = [f FF(2:end)];
H = [H1; H(2:end,:)];
Amp = abs(H);
Phase = 180/pi*unwrap(angle(H));
figure()
subplot(2,1,1)
plot(f,20*log10(Amp),'--o')
xlabel('f [Hz]')
ylabel('Gain mm/N [dB]')
subplot(2,1,2)
plot(f,Phase,'--o')
xlabel('f [Hz]')
ylabel('Phase [deg]')
legend('m_1','m_2','m_3','m_4')

zeta = [0.0126, 0.0076, 0.004, 0.002];
for ii = 1:4
    load(['bode_FIRST_ZERO(13Hz)Res' num2str(ii) '.mat'])
    figure()
    subplot(2,1,1)
    plot(FF,20*log10(abs(H)),'--o')
    xlabel('f [Hz]')
    ylabel('Gain mm/N [dB]')
    subplot(2,1,2)
    plot(FF,180/pi*unwrap(angle(H)),'--o')
    xlabel('f [Hz]')
    ylabel('Phase [deg]')
    legend('m_1','m_2','m_3','m_4','location','best')
    [~,res] = max(abs(H(:,1)));
    disp(['F_d = ' num2str(FF(res)) ' F_n = ' num2str(FF(res)/sqrt(1-zeta(ii)^2))]);

    [mH, mind] = max(abs(H));
    [mmH, mmind] = max(mH);
    Norm_Max_H = H(mind(mmind),:)./abs(H(mind(mmind),mmind));
    v = sign(imag(Norm_Max_H)).'.*abs(Norm_Max_H).'
    figure()
    polar(angle(H),abs(H),'-o')
    hold on
    for jj = 1:4
        polar(angle(H(mind(jj),jj)),abs(H(mind(jj),jj)),'*')
    end
    legend('m_1','m_2','m_3','m_4','location','best')
    hold off
end