%  Matlab code to estimate the settle time for a sine input
%  Shai  16-Oct-2015
clc;zeta = 0.03; wn = 2*pi*31;  % parameters of oscillator to simulate)
[A,B,C,D] = ord2(wn,zeta); 	% create a state-space model of m,k,c
SDOF = ss(A,B,C,D);         % convert to a state space object
T =3;                    % total time
U = 1;                      % Sine amplitude
freq_vec=5:5:100;
Fs = 500;           % sampling rate Hz
N = Fs*T;           % No. of collected data points
t =(0:N)*1/Fs;    % time vector
Estimation_Error=0;t_SS_Avg=0;
for f_excitation=freq_vec
    w = 2*pi*f_excitation;
    u = U*sin(w*t);     % create sine input
    y = lsim(SDOF,u,t);    % run the simulation
    
    % Find the envelop efficiently using 'reshape' and 'max' commands
    Nsam=round(1/f_excitation*Fs); % # of samples per cycle in the response
    % reshape y to find the max in each period
    y_cut=y(1:end-mod(length(y),Nsam));
    tmp=reshape(y_cut,Nsam,length(y_cut)/Nsam);
    Envelop=max(tmp);
    t_env=linspace(0,max(t),length(Envelop));
    % Fint the time to steady state
    idx_SS=find(abs(Envelop-Envelop(end))/Envelop(end)<1e-4,1);
    t_SS=t_env(idx_SS);
    t_est=1/zeta/wn/sqrt(1-zeta^2)*log(1e4);
    
    plot(t,y,t_env,Envelop,'k--','linewidth',4);set(gca,'fontsize',14);
    title({['Response to sin @ ',num2str(f_excitation),'Hz input'],...
        ['Time to Steady-state is ',num2str(t_SS,'%.2f'),'s'],...
        ['Estimated settle time is ',num2str(t_est,'%.2f'),'s']},'fontsize',14);
    xlabel 'Time (s)'
    ylabel 'Response (a.u.)'
    drawnow;
    % Display the difference btw estimation and reality
    Estimation_Error=Estimation_Error+...
        round(abs(t_SS-t_est)/t_SS*100)/length(freq_vec);
    t_SS_Avg=t_SS_Avg+t_SS/length(freq_vec);
    pause(.5);
end
% Display the estimation error
Estimation_Error
% Average settle time
t_SS_Avg