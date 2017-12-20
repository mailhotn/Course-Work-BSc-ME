

Fs=500; %sampling frequency
t=0:1/Fs:1; %time base - upto 1 second

f0=1;% starting frequency of the chirp
f1=Fs/20; %frequency of the chirp at t1=1 second

t0=t(1);
T=1-t0;
k=(f1-f0)/T;
x=cos(2*pi*(k/2*t+f0).*t); % same as x = chirp(t,f0,1,f1);

subplot(2,2,1);
plot(t,x,'k');title('Chirp Signal');xlabel('Time(s)');ylabel('Amplitude');

% Compute the analytic signal and differentiate its phase to measure the instantaneous frequency. The scaled derivative yields a meaningful estimate
subplot(2,2,2);
z=hilbert(x);
instfreq = Fs/(2*pi)*diff(unwrap(angle(z)));
plot(t(2:end),instfreq,t(2:end),(k*t(2:end)+f0));
title('Hilbert transform');xlabel('Time(s)');ylabel('inst freq');

% fft
subplot(2,2,3);
L=length(x);
NFFT = 2^nextpow2(L);
Y = fft(x,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
plot(f,2*abs(Y(1:NFFT/2+1))) 
xlim([0 50]);
title('FFT');xlabel('Frequency (Hz)');ylabel('Amplitude');
