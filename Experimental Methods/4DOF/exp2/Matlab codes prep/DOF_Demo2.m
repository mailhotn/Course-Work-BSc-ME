clc;clf;disp('Press F10 to advance line by line');

button = questdlg('Clear all variables?');
switch button
    case 'Yes'
        clear all;
    case 'Cancel'
end

dbstop at 23;

str={'Chirp example',...
    'Chirp example 2',...
    'Hilbert example',...
    'Stepped sine example',...
    'Fit frf example'};

[s,v] = listdlg('PromptString','Select a file:',...
    'SelectionMode','single',...
    'ListString',str);
if ~v,return;end;

switch str{s}
    case 'Chirp example'
        % Chirp Example.m
        Fs = 5000;          % sampling rate Hz
        T = 20;             % total time
        N = Fs*T;           % No. of collected data points
        t = [0:N]'*1/Fs;    % time vector
        
        f0 = 10; f1 = 40;   % start and end frequencies
        
        phi = 2*pi*(f0*t+(f1-f0)/2/T*t.^2); % phase yielding d(phi)/dt=f0+t/T*(f1-f0)
        f = f0 +(f1-f0)/T*t;
        u = sin(phi);                   % create chirp (Swept sine), unit amplitude
        zeta = 0.03; wn = 2*pi*30;     	% parameters of oscillator to simulate)
        [A,B,C,D] = ord2(wn,zeta);      % create a state-space model of m,k,c
        SDOF = ss(A,B,C,D);             % convert to a state space object
        y = lsim(SDOF,u,t);             % run the simulation
        %I_f = gradient(phi)./gradient(t)/2/pi; % instantaneous frequency
        t_wn = T*(wn/2/pi-f0)/(f1-f0);
        
        figure(1);
        subplot(3,1,[1 2])
        [AX, ~, H2] = plotyy(t,y,t,f);  % plot the input and output
        hold on
        set(H2,'LineWidth',2);
        xlim(AX(1),[t(1) t(end)])
        xlim(AX(2),[t(1) t(end)])
        % xlabel('$time$[sec]','interpreter','latex','FontSize',18)
        ylabel(AX(1),'$y(t)$','interpreter','latex','FontSize',18)
        ylabel(AX(2),'$instantaneous$ $frequency$[Hz]','interpreter','latex','FontSize',18)
        title('$System$ $Response$','interpreter','latex','FontSize',18)
        set(AX(1),'FontSize',14);
        set(AX(2),'FontSize',14);
        set(gcf,'Color','w');
        
        t_wn = 1/((f1-f0)/T/(wn/2/pi-f0));
        
        subplot(3,1,3)
        h = plot(t,u,[t_wn t_wn],[-1 1]);
        set(h(2),'LineWidth',4,'Color','r')
        xlabel('$time$[sec]','interpreter','latex','FontSize',18)
        ylabel('$u(t)$','interpreter','latex','FontSize',18)
        title('$Input$ $force$','interpreter','latex','FontSize',18)
        set(gca,'FontSize',14);
        
    case 'Chirp example 2'
        
        TT = [1 3 10 66  500]; % total time
        TN = length(TT);
        HL = zeros(TN,1);
        for kk = 1:TN
            
            T = TT(kk);         % total time
            Fs = 500;           % sampling rate Hz
            
            N = Fs*T;           % No. of collected data points
            t = [0:N]'*1/Fs;    % time vector
            
            f0 = 24; f1 = 34;   % start and end frequencies
            phi = 2*pi*(f0*t+(f1-f0)/2/T*t.^2);
            f = f0 +(f1-f0)/T*t;
            u = sin(phi);             	% create chirp (Swept sine)
            zeta = 0.03; wn = 2*pi*30; 	% parameters of oscillator to simulate)
            [A,B,C,D] = ord2(wn,zeta); 	% create a state-space model of m,k,c
            SDOF = ss(A,B,C,D);       	% convert to a state space object
            y = lsim(SDOF,u,t);         % run the simulation
            I_f = gradient(phi)./gradient(t)/2/pi; % instantaneous frequency
            
            figure(1)
            subplot(TN,1,kk)
            plot(t ,y)
            [AX, ~, H2] = plotyy(t,y,t,I_f);  % plot the output
            hold on
            set(H2,'LineWidth',2);
            xlim(AX(1),[t(1) t(end)])
            xlim(AX(2),[t(1) t(end)])
            xlabel('time[sec]','interpreter','latex','FontSize',12)
            ylabel(AX(1),'$y(t)$','interpreter','latex','FontSize',12)
            ylabel(AX(2),'$f(t)$[Hz]','interpreter','latex','FontSize',12)
            title(sprintf('$T = %d$',T),'interpreter','latex','FontSize',12)
            set(AX(1),'FontSize',10);
            set(AX(2),'FontSize',10);
            set(gcf,'Color','w');
            hold off
            
            figure(2);subplot(TN,1,kk);
            HL(kk) = plot(f,abs(hilbert(y)),'LineWidth',2);
            hold all
            xlabel('$f$[Hz]','interpreter','latex','FontSize',18)
            ylabel('$Instantaneous amplitude$','interpreter','latex','FontSize',18)
            set(gcf,'Color','w');
            set(gca,'FontSize',18);
            
            figure(3);subplot(TN,1,kk);
            %             http://www.mathworks.com/help/signal/ug/hilbert-transform-and-instantaneous-frequency.html
            z=hilbert(y);
            instfreq = Fs/(2*pi)*diff(unwrap(angle(z)));
            plot(f(2:end),abs(log(instfreq-f(2:end))));
            title('Instantaneous frequency');
            xlabel('Frequency (Hz)');ylabel('Response');
            
            figure(4);subplot(TN,1,kk);
            L=length(y);
            NFFT = 2^nextpow2(L);
            Y = fft(y,NFFT)/L;
            f = Fs/2*linspace(0,1,NFFT/2+1);
            plot(f,2*abs(Y(1:NFFT/2+1)));
            xlabel('Frequency (Hz)');ylabel('Response');
            title('FFT');
            
        end
        figure(2)
        plot([wn/2/pi wn/2/pi],[min(abs(hilbert(y))) max(abs(hilbert(y)))],'--r','LineWidth',2);
        hold off
        hl = legend(HL,sprintf('T=%d',TT(1)),sprintf('T=%d',TT(2)),sprintf('T=%d',TT(3)),sprintf('T=%d',TT(4))...
            ,sprintf('T=%d',TT(5)));
        set(hl,'interpreter','latex','FontSize',14);
        
    case 'Hilbert example'
        
        Fs=500; %sampling frequency
        t=0:1/Fs:1; %time base - upto 5 second
        
        f0=10;% starting frequency of the chirp
        f1=Fs/20; %frequency of the chirp at t1=1 second
        
        t0=t(1);
        T=1-t0;
        k=(f1-f0)/T;
        
        %         Example #1 same as x = chirp(t,f0,1,f1);
        phi=2*pi*(k/2*t+f0).*t;x=cos(phi); % Example #1
        
        %         Example #2 exponential chirp
        %         phi=2*pi*(k/2*(1-exp(-t))+f0).*t;x=exp(-(t-t(end)/2).^2).*cos(phi);
        
        %         Example #3 perfect cosine
        %                 phi=2*pi*f0*t;x=cos(phi);
        
        %         Example #4 two cosines
        %         phi=2*pi*f0*t;x=cos(phi)+.5*cos(1.732*phi);
        
        %         Example #5 cosine with noise
        %         phi=2*pi*f0*t;x=cos(phi)+rand(1,length(t));
        
        subplot(2,2,1);
        plot(t,x,'k');title('Signal');xlabel('Time(s)');ylabel('Amplitude');
        ylim([-1.5 1.5]);
        
        % Compute the analytic signal and differentiate its phase to measure the instantaneous frequency. The scaled derivative yields a meaningful estimate
        subplot(2,2,2);
        z=hilbert(x);
        instfreq = Fs/(2*pi)*diff(unwrap(angle(z)));
        plot(t(2:end),instfreq,t(2:end),diff(phi)/2/pi*Fs,'linewidth',4);
        legend('angle(z)','dphi/dt');
        title('Hilbert transform inst. freq');xlabel('Time(s)');ylabel('inst freq');
        
        subplot(2,2,3);
        plot(t,x,t,abs(z),'linewidth',4);
        legend('Signal','abs(z)');ylim([-1.5 1.5]);
        title('Hilbert transform inst amplitude');xlabel('Time(s)');ylabel('inst amplitude');
        
        % fft
        subplot(2,2,4);
        L=length(x);
        NFFT = 2^nextpow2(L);
        Y = fft(x,NFFT)/L;
        f = Fs/2*linspace(0,1,NFFT/2+1);
        plot(f,2*abs(Y(1:NFFT/2+1)))
        xlim([0 50]);
        title('FFT');xlabel('Frequency (Hz)');ylabel('Amplitude');
        
        
    case  'Stepped sine example'
        
        % Stepped sine - 2 points
        % clear  G PS
        p = 1;
        f = linspace(1,40,100);     % frequencies vector
        zeta = 0.03; wn = 2*pi*30;  % parameters of oscillator to simulate)
        [A,B,C,D] = ord2(wn,zeta); 	% create a state-space model of m,k,c
        SDOF = ss(A,B,C,D);         % convert to a state space object
        T = 100;                    % total time
        U = 1;                      % Sine amplitude
        
        % Stepped Sine base on 2 points - No noise
        noise = 0;  % when 0, no measurement noise, set noise =1e-4
        for w = 2*pi*f
            
            Fs = 500;           % sampling rate Hz
            % total time
            N = Fs*T;           % No. of collected data points
            t = [0:N]'*1/Fs;    % time vector
            
            u = U*sin(w*t);     % create sine input
            
            y = lsim(SDOF,u,t)+noise*randn(size(t));    % run the simulation
            i1 = find(t>0.8*max(t));                    % wait 80% of total time hopping we have steady state
            i1 = [N N+1]; %<<<<<<<<<- remove comment to enforce 2 point only based solution
            t1 =t(i1);  y1 = y(i1); u1 = u(i1);         % all data after 80% of time passed
            cs = [cos(w*t1) sin(w*t1)];                 % form regression matrix
            CD = cs\y1;                                 % compute the least squares solution of C,D
            AB = cs\u1;                                 % compute the least square of A,B foc A*cos(w*t)+B*sin(w*t)
            
            figure(1)
            G(p,1) = norm(CD)/norm(AB);                 % compute gain (sqrt(C^2+D^2)/sqrt(A^2+B^2))
            PS(p,1) = atan2(CD(1),CD(2))-atan2(AB(1),AB(2)); % compute phase
            subplot(2,1,1)
            plot(f(1:p),G(1:p),'o','MarkerSize',6)      % plot gain for all frequencies until point p
            xlabel('$f$[Hz]','interpreter','latex','FontSize',18)
            ylabel('Gain','interpreter','latex','FontSize',18)
%             title('2 Points - No Noise','interpreter','latex','FontSize',18)
            set(gca,'FontSize',14)
            set(gcf,'Color','w')
            
            subplot(2,1,2)
            plot(f(1:p),180/pi*PS(1:p),'o') % plot phase in degress for all frequencies until point p
            ylabel('Phase[Deg]','interpreter','latex','FontSize',18)
            xlabel('$f$[Hz]','interpreter','latex','FontSize',18)
            set(gca,'FontSize',14)
            drawnow   % force refresh of display
            
            p = p+1;  % next frequency from vector
            shg
        end
        
        % Stepped Sine base on 2 points - with noise
        p = 1;
        noise = 1e-4;  % when 0, no measurement noise, set noise =1e-4
        for w = 2*pi*f
            
            Fs = 500;           % sampling rate Hz
            % total time
            N = Fs*T;           % No. of collected data points
            t = [0:N]'*1/Fs;    % time vector
            
            u = U*sin(w*t);     % create sine input
            
            y = lsim(SDOF,u,t)+noise*randn(size(t));    % run the simulation
            i1 = find(t>0.8*max(t));                    % wait 80% of total time hopping we have steady state
            i1 = [N N+1]; %<<<<<<<<<- remove comment to enforce 2 point only based solution
            t1 = t(i1);  y1 = y(i1); u1 = u(i1);         % all data after 80% of time passed
            cs = [cos(w*t1) sin(w*t1)];                 % form regression matrix
            CD = cs\y1;                                 % compute the least squares solution of C,D
            AB = cs\u1;                                 % compute the least square of A,B foc A*cos(w*t)+B*sin(w*t)
            
            figure(2)
            G(p,1) = norm(CD)/norm(AB);                 % compute gain (sqrt(C^2+D^2)/sqrt(A^2+B^2))
            PS(p,1) = atan2(CD(1),CD(2))-atan2(AB(1),AB(2)); % compute phase
            subplot(2,1,1)
            plot(f(1:p),G(1:p),'o','MarkerSize',6)      % plot gain for all frequencies until point p
            xlabel('$f$[Hz]','interpreter','latex','FontSize',18)
            ylabel('Gain','interpreter','latex','FontSize',18)
%             title('2 points - With Noise','interpreter','latex','FontSize',18)
            set(gca,'FontSize',14)
            set(gcf,'Color','w')
            
            subplot(2,1,2)
            plot(f(1:p),180/pi*PS(1:p),'o') % plot phase in degress for all frequencies until point p
            ylabel('Phase[Deg]','interpreter','latex','FontSize',18)
            xlabel('$f$[Hz]','interpreter','latex','FontSize',18)
            set(gca,'FontSize',14)
            drawnow   % force refresh of display
            
            p = p+1;  % next frequency from vector
            shg
        end
        
        % Stepped Sine N points - with noise
        p = 1;
        noise = 1e-4;  % when 0, no measurement noise, set noise =1e-4
        for w = 2*pi*f
            
            Fs = 500;           % sampling rate Hz
            % total time
            N = Fs*T;           % No. of collected data points
            t = [0:N]'*1/Fs;    % time vector
            
            u = U*sin(w*t);     % create sine input
            
            y = lsim(SDOF,u,t)+noise*randn(size(t));    % run the simulation
            i1 = find(t>0.8*max(t));                    % wait 80% of total time hopping we have steady state
            %     i1 = [N N+1]; %<<<<<<<<<- remove comment to enforce 2 point only based solution
            t1 = t(i1);  y1 = y(i1); u1 = u(i1);         % all data after 80% of time passed
            cs = [cos(w*t1) sin(w*t1)];                 % form regression matrix
            CD = cs\y1;                                 % compute the least squares solution of C,D
            AB = cs\u1;                                 % compute the least square of A,B foc A*cos(w*t)+B*sin(w*t)
            
            figure(3)
            G(p,1) = norm(CD)/norm(AB);                 % compute gain (sqrt(C^2+D^2)/sqrt(A^2+B^2))
            PS(p,1) = atan2(CD(1),CD(2))-atan2(AB(1),AB(2)); % compute phase
            subplot(2,1,1)
            plot(f(1:p),G(1:p),'o','MarkerSize',6)      % plot gain for all frequencies until point p
            xlabel('$f$[Hz]','interpreter','latex','FontSize',18)
            ylabel('Gain','interpreter','latex','FontSize',18)
%             title('N points - With Noise','interpreter','latex','FontSize',18)
            set(gca,'FontSize',14)
            set(gcf,'Color','w')
            
            subplot(2,1,2)
            plot(f(1:p),180/pi*PS(1:p),'o') % plot phase in degress for all frequencies until point p
            ylabel('Phase[Deg]','interpreter','latex','FontSize',18)
            xlabel('$f$[Hz]','interpreter','latex','FontSize',18)
            set(gca,'FontSize',14)
            drawnow   % force refresh of display
            
            p = p+1;  % next frequency from vector
            shg
        end
        
    case 'Fit frf example'
        
        % fit_frf_danciger
        clc; clear all; close all;
        % simulate experiment 2 DOF
        % define masses and springs' stiffness
        m1 = 1; m2 = 1;
        k1 = 10; k2 =3 ;
        % define mass, damping and stiffness marices
        M = diag([m1 m2]); K = [2*k1+k2 -k2; -k2 2*k1+k2];
        C = [.1 0;0 .1];
        % define State-Space marices
        A = [zeros(2) eye(2); -M\[K C]];  B = [zeros(2);inv(M)];
        c = [eye(2) zeros(2)]; D=zeros(2);
        
        sys = ss(A,B,c,D);
        dt = 1/100;             % time step
        N = 25000;              % No. of collected data points
        Fs = 1/dt;              % Sampling frequency
        t = [0:N-1]'*dt;        % time vector
        [wn1,zn] = damp(sys);   % get damping and wn
        df=Fs/N;
        
        [vv,w2] = eig(K,M);
        wn = sqrt(diag(w2));    % find natural frequencies
        % stepped sine
        
        w = linspace(wn(2)*0.95,wn(2)*1.05,50);
        i1 = find(t>5/(zn(end)*wn1(end)));      % choose times where we've reached steady-state
        i2= find(t>t(i1(1)) & t<t(i1(1))+5);    % range for plotting
        
        p=1;
        for w1=w
            u = sin(w1*t);                      % excitation
            y = lsim(sys(:,1),u,t);
            
            figure(1)
            plot(t,y)
            xlabel('time[sec]','interpreter','latex','FontSize',18)
            ylabel('Output','interpreter','latex','FontSize',18)
            set(gcf,'Color','w')
            set(gca,'FontSize',14)
            drawnow, shg
            pause(.1) % wait 0.1 sec
            
            %fit sine
            cs = [cos(w1*t) sin(w1*t) ones(size(t))]; %y=ac*cos(wt)+as*sin(wt)+dc
            au = cs(i1,:)\u(i1,:);      % fit input cos and sine coefficients
            ay = cs(i1,:)\y(i1,:);      % fit out cos/sin
            
            figure(3)
            plot(t(i2),y(i2,:),'.',t(i2),cs(i2,:)*ay)
            xlabel('time[sec]','interpreter','latex','FontSize',18)
            ylabel('Output','interpreter','latex','FontSize',18)
            hl = legend('$y_1$','$y_2$','fitted-$y_1$','fitted-$y_2$');
            set(hl,'interpreter','latex','FontSize',14)
            title(sprintf('$f =%3.3f$ [Hz]',w1/2/pi),'interpreter','latex','FontSize',14)
            set(gcf,'Color','w')
            set(gca,'FontSize',14)
            pause(.1), drawnow, shg
            
            H(p,:) = (ay(1,:)-1i*ay(2,:))/(au(1,:)-1i*au(2,:)); % compute H(i*w1)
            
            figure(2)
            plot(w(1:p),abs(H(1:p,:)),'o-','LineWidth',2,'MarkerSize',6)
            xlabel('$\omega$[Rad/s]','interpreter','latex','FontSize',18)
            ylabel('$|H(\omega)|$','interpreter','latex','FontSize',18)
            set(gcf,'Color','w')
            set(gca,'FontSize',14)
            drawnow, shg
            
            figure(4)
            plot(u(i2),y(i2,:),'LineWidth',2)
            xlabel('force','interpreter','latex','FontSize',18)
            ylabel('response','interpreter','latex','FontSize',18)
            hl = legend('$y_1$','$y_2$');
            set(hl,'interpreter','latex','FontSize',14)
            title(sprintf('phases: $H_{1,1}=%3.3f$[deg],  $H_{2,1}=%3.3f$ [deg]',angle(H(p,:))*180/pi),'interpreter','latex','FontSize',14)
            axis equal
            set(gcf,'Color','w')
            set(gca,'FontSize',14)
            drawnow, shg
            pause(.2)
            p=p+1;
            tilefigs
        end
        
        % plot the frequency response as a Nyquist plot
        figure(5)
        plot(H,'LineWidth',2)
        xlabel('Real${H}$','interpreter','latex','FontSize',18)
        ylabel('Imaginary${H}$','interpreter','latex','FontSize',18)
        axis equal
        set(gcf,'Color','w')
        set(gca,'FontSize',14);
        
end