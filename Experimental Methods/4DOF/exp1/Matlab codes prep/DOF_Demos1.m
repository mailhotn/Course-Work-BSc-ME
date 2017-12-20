clc;clf;disp('Press F10 to advance line by line');

button = questdlg('Clear all variables?');
switch button
    case 'Yes'
        clear all;
    case 'Cancel'
end

dbstop at 28;

str={'Quantization demo',...
    'Fit graph using fit command',...
    'Fourier demo',...
    'Multi DOF example',...
    'Sine decay example',...
    'Show fourier of measured data',...
    'Multi DOF eigenvalue example',...
    'Line fit demo',...
    'Sine fit demo',...
    'Stroboscope demo'};

[s,v] = listdlg('PromptString','Select a file:',...
    'SelectionMode','single',...
    'ListString',str);
if ~v,return;end;

switch str{s}
    case 'Quantization demo'
        % Quantization DEMO
        clear all; close all; clc
        
        Vs = 5; N = 5;
        
        Nc = 15.67; Np = 1021;
        wt = linspace(0,2*pi*Nc,Np)';
        y = Vs*sin(wt);
        
        for N = 3:5
            
            Q = 2*Vs/2^N;
            yq = quant(y,Q);
            
            figure(N)
            subplot(2,1,1)
            plot(wt,y,'LineStyle','none','Marker','o','MarkerSize',5,'MarkerFaceColor','r','MarkerEdgeColor','r');
            hold on
            plot(wt,yq,'LineStyle','none','Marker','o','MarkerSize',5,'MarkerFaceColor','b','MarkerEdgeColor','b');
            lh = legend('Unquantized','Quantized');
            set(lh,'interpreter','latex','FontSize',18)
            title(sprintf('Bits = %d',N),'interpreter','latex','FontSize',18)
            set(gcf,'Color','w');
            set(gca,'FontSize',14);
            
            subplot(2,1,2)
            plot(wt,y,'LineStyle','none','Marker','o','MarkerSize',5,'MarkerFaceColor','r','MarkerEdgeColor','r');
            hold on
            plot(wt,yq,'LineStyle','none','Marker','o','MarkerSize',5,'MarkerFaceColor','b','MarkerEdgeColor','b');
            ylim([-(N-1)*Q,(N-1)*Q])
            axh = gca;
            set(axh,'YTick',[-(N-1)*Q:Q:(N-1)*Q])
            set(axh,'XTick',0)
            lh = legend('Unuantized','Quantized');
            set(lh,'interpreter','latex','FontSize',18)
            grid on
            set(gca,'FontSize',14);
        end
        
    case 'Fit graph using fit command'
        
        wn = 2*pi*12;  zeta = 0.01; sigma = zeta*wn; % true parameters
        Fs = 200;                                    % 200 samples/sec
        N = 1000; dt = 1/Fs;
        t = (0:N-1)'*dt;
        
        [A,B,C,D]=ord2(wn,zeta);  sys = ss(A,B,C,D);    % 1dof oscillator model
        x0 = 1; v0 = 0.2;                            	% Initial conditions
        y = initial(sys,[x0;v0],t);                     % simulate experiment
        
        % add noise
        y = y+0.05*randn(size(y));
        
        % now fit a model using Matlab's curve fitting toolbox
        ft = fittype('a+b*exp(c*x)*cos(d*x+e)'); % Asuumed model
        %coeffnames(ft)
        
        [crv,gof,output] = fit(t,y,ft,'StartPoint',[0,1/2,-.2,70,0.2]);
        yf = feval(crv,t);
        
        figure(1);
        h = plot(t,[y yf output.residuals]);
        set(h,'LineWidth',2)
        set(h(1),'Marker','x')
        set(h(3),'Color','r')
        xlabel('time[sec]','interpreter','latex','FontSize',18)
        ylabel('$y$','interpreter','latex','FontSize',18)
        shg
        hl = legend('measured','fitted','error');
        set(hl,'interpreter','latex','FontSize',18)
        set(gcf,'Color','w')
        set(gca,'FontSize',14)
        
        title(sprintf('$f_n=%.1f$[Hz], $\\zeta=%.4f$',crv.d/2/pi,-crv.c/crv.d),'interpreter','latex','FontSize',18)
        
    case 'Fourier demo'
        
        dt = 1/5000;        % time step
        N = 5000;           % Number of samples
        Fs = 1/dt;          % Sampling frequency
        t = [0:N-1]'*dt;    % time vector
        
        df = Fs/N;          % Frequency resolution
        
        f1 = 20*df;
        f2 = 65*df;
        
        x = 2*sin(2*pi*f1*t)+3*cos(2*pi*f2*t);  % The input signal
        
        figure(1)
        stem([f1;f2],[2;3],'LineWidth',2,'MarkerFaceColor',lines(1))
        xlim([0,80])
        xlabel('$f$[Hz]','interpreter','latex','FontSize',18)
        ylabel('$|X(\omega)|$','interpreter','latex','FontSize',18)
        title('The frequncies comprising the input signal','interpreter','latex','FontSize',18)
        set(gca,'FontSize',14)
        set(gcf,'Color','w')
        shg
        
        figure(2)
        [X,f] = FourierSeriesDFT(x,Fs);
        plot(f,abs(X),'o-','LineWidth',2,'MarkerFaceColor',lines(1))
        xlim([0,80])
        xlabel('$f$[Hz]','interpreter','latex','FontSize',18)
        ylabel('$|X(\omega)|$','interpreter','latex','FontSize',18)
        title('DFT','interpreter','latex','FontSize',18)
        set(gca,'FontSize',14)
        
    case 'Multi DOF example'
        
        % mdof_matlab_example
        close all; clear all; clc;
        
        m1 = 1; m2 = 1;     % Masses
        k1 = 10; k2 = 3;    % Springs
        
        M = diag([m1 m2]);              % Mass matrix
        K = [2*k1+k2 -k2; -k2 2*k1+k2]; % Stiffness matrix
        C = [.1 0;0 .1];                % Damping matrix
        
        % State-Space
        A = [zeros(2) eye(2); -M\[K C]];  B = [zeros(2);inv(M)];
        c =[eye(2) zeros(2)]; D=zeros(2);
        
        sys = ss(A,B,c,D);
        
        step(sys)   % compute the step response
        set(gcf,'color','w')
        pause
        
        bode(sys)   % compute bode
        set(gcf,'color','w')
        
        % compute eigenvectors and natural frequencies
        [V,lam]=eig(K,M);
        
        wn = sqrt(diag(lam));
        
        w=wn(1);
        x1 =(K-w^2*M+1i*w*C)\[1;0];
        x1/x1(1)
        
        w=wn(2);
        x2 =(K-w^2*M+1i*w*C)\[1;0];
        x2/x2(1)
        
    case 'Sine decay example'
        
        ff=[9.6  19.3]; zz=[0.01 0.03]
        
        for q=1:2
            
            wn = 2*pi*ff(q);  zeta = zz(q); sigma = zeta*wn; % true parameters
            Fs = 200;                                    % 200 samples/sec
            N = 1000; dt = 1/Fs;
            time = (0:N-1)'*dt;
            
            [A,B,C,D]=ord2(wn,zeta);  sys = ss(A,B,C,D);    % 1dof oscillator model
            x0 = 1; v0 = 0.2;                            	% Initial conditions
            y = initial(sys,[x0;v0],time);                     % simulate experiment
            
            % add noise
            data = y+0.01*randn(size(y));
            
            plot(time,data)
            title('choose 2 points on exp. decay part')
            shg
            [xx yy]=ginput(2);
            
            xx=sort(xx);
            i1 = find(time>=xx(1) & time<=xx(2));
            
            plot(time(i1),data(i1,:))
            t=time(i1);
            
            y=detrend(data);
            y3= y(i1);
            plot(time,y,time(i1),y(i1))
            % fit decaying sine
            drawnow
            
            ft = fittype('a+b*exp(c*x)*cos(d*x+e)');
            %coeffnames(ft)
            
            [crv,gof,output] = fit(t,y3,ft,'StartPoint',[0,1/2,-.2,ff(q)*2*pi,0.2]);
            
            yf = feval(crv,t);
            
            plot(t,y3,t,yf ,t,output.residuals)
            legend measured fitted error
            
            title(sprintf('f_d=%.1f Hz, f_n=%.1f Hz, zeta=%.4f',crv.d/2/pi, sqrt(crv.c^2+crv.d^2)/2/pi,-crv.c/sqrt(crv.d^2+crv.c^2)))
            Fd = crv.d/2/pi;
            FN = sqrt(crv.c^2+crv.d^2)/2/pi;
            ZETA = -crv.c/sqrt(crv.d^2+crv.c^2);
        end
        
    case 'Show fourier of measured data'
        
        % load saved measurement
        SensorNo = 1:4;
        
        cc = hsv(4);
        h1 = zeros(4,1);
        for q = 3 %1:4
            
            load('D3.mat')
            % choose data range  to be fitted
            dt = t(2)-t(1); Fs=1/dt;
            
            for kk = SensorNo
                y = ch(:,kk);
                [Y,f] = FourierSeriesDFT(y.*hanning(length(y)),Fs);
                
                Ydb = 20*log10(abs(Y));
                
                pks = find(Ydb>-66);
                
                figure
                h = plot(f,Ydb,'--o',f(pks),Ydb(pks),'x');
                set(h(1),'Color',lines(1),'MarkerFaceColor',lines(1),'MarkerEdgeColor',lines(1))
                set(h(2),'MarkerFaceColor','r','MarkerSize',10,'LineWidth',3)
                xlabel('$f$[Hz]','Interpreter','Latex','FontSize',18)
                ylabel('$|Y(\omega)|$[dB]','Interpreter','Latex','FontSize',18)
                title(sprintf('Measurements of mass %d',q),'Interpreter','Latex','FontSize',18)
                xlim([2 80])
                set(gca,'FontSize',14)
                set(gcf,'Color','w')
                shg
                title(sprintf('DFT measurments of mass No. %d',kk),'Interpreter','Latex','FontSize',18);
            end
        end
        
    case 'Multi DOF eigenvalue example'
        % mdof_matlab_example
        
        m1 = 1; m2 = 1;     % Masses
        k1 = 10; k2 = 3;    % Springs
        
        M = diag([m1 m2]);              % Mass matrix
        K = [2*k1+k2 -k2; -k2 2*k1+k2]; % Stiffness matrix
        C = [.1 0;0 .1];                % Damping matrix
        
        % State-Space
        A = [zeros(2) eye(2); -M\[K C]];  B = [zeros(2);inv(M)];
        c =[eye(2) zeros(2)]; D=zeros(2);
        
        sys = ss(A,B,c,D);
        
        % eigenvalues
        
        [V,W2] = eig(K,M);      % solve generalized eigenvalue problem
        
        wn = sqrt(diag(W2));    % get natural frequencies
        two_zeta_wn = V.'*C*V;  % true for proportinal damping
        zeta = diag(two_zeta_wn)./wn/2;
        
        % show orthogonality
        mr = V.'*M*V;  % should be unity matrix
        kr = V.'*K*V;  % should be  diagonal with wn^2 as entries
        
        N = size(V,1); nw = 1000;
        w = linspace(0,max(wn)*2,nw)';
        H = zeros(nw,N,N); Hmode = zeros(nw,N,N,N);
        for n1 = 1:nw,
            for n2 = 1:N
                den = wn(n2)^2-w(n1).^2+1i*two_zeta_wn(n2,n2)*w(n1);
                H(n1,:,:)=squeeze(H(n1,:,:))+(V(:,n2)*V(:,n2).')/den;
                Hmode(n1,:,:,n2)=squeeze(Hmode(n1,:,:,n2))+(V(:,n2)*V(:,n2).')/den;
            end
        end
        
        [Mag]=bode(sys(1,1),w);
        figure
        h = plot(w/2/pi,squeeze(Mag),'o',w/2/pi,abs(squeeze(H(:,1,1))),...
            w/2/pi,abs(squeeze(Hmode(:,1,1,1))),'-',...
            w/2/pi,abs(squeeze(Hmode(:,1,1,2))),'-');
        set(h,'LineWidth',2)
        set(h(1),'MarkerSize',4)
        hl = legend(h,'$bode$','$mode$ $sum$','$mode_1$','$mode_2$');
        set(hl,'interpreter','latex','FontSize',18)
        xlim([w(1)/2/pi w(end)/2/pi])
        shg
        set(gcf,'color','w')
        set(gca,'FontSize',14)
        xlabel('frequency','interpreter','latex','FontSize',18)
        ylabel('$|H_{1,1}(\omega)|$','interpreter','latex','FontSize',18)
        
    case 'Line fit demo'
        % line_fit_demo.m
        
        x = [1 2 3 4 6 7]';                 % measured x values
        y = [3 5.1 7.3 8.9 12.9 15.03]';    % corresponding y values
        
        
        % assuming y=y(x) and y=a*x+b +e   (e - error)
        % we can show that e= y-(a*x+b)
        % therefore the sum of squares is e'*e
        %
        % rewrite the model as:
        %   [y(1)]  [x(1) 1]
        %   [y(2)]  [x(2) 1]
        %   [y(3)] =[x(3) 1] [a]
        %   [y(4)] =[x(4) 1] [b]
        %
        % or
        % y = A*q+ B,          q=[a;b]
        % e = error = y-(A*q+B),
        % minimizing e'*e -->    (A'*A)*q = A'*y -> q=inv(A'*A)*A'*y
        
        % solution method #1
        A = [x ones(size(x))];
        q = inv(A'*A)*A'*y;
        e = y-A*q;   ErrorFun=e'*e;
        
        % solution method #2 - using Matlab's bulitin function
        q1 = polyfit(x,y,1);
        
        % solution method #3 - unsing the curve fitting tollbox
        [q2 o] = fit(x,y,'a*x+b')
        
        % since all the methods produced the same result, we plot one
        
        yfit1 = A*q;  yfit2=polyval(q1,x);
        figure(1);
        h = plot(x,[y yfit1 y-yfit1]);
        set(h(1),'Marker','o','MarkerSize',6,'LineStyle','none')
        set(h(2),'LineWidth',2)
        set(h(3),'Marker','o','MarkerSize',6,'LineStyle','none','Color','r')
        xlabel('$x$','interpreter','latex','FontSize',18)
        ylabel('$y$','interpreter','latex','FontSize',18)
        hl = legend('$measured$','$fitted$','Error');
        set(hl,'interpreter','latex','FontSize',18)
        set(gcf,'Color','w')
        set(gca,'FontSize',14)
        figure(gcf)
        
    case 'Sine fit demo'
        
        % fit_sine_danciger
        % simulate experiment
        dt = 1/600;             % time step
        N = 250;                % No. of collected data points
        Fs = 1/dt;              % Sampling frequency
        t = [0:N-1]'*dt;        % time vector
        w = 2*pi*20.123;        % arbitrary frequency Rad/s
        phi = 22*pi/180;        % arbitrary phase Rad
        A = 1.2345;             % amplitude   V
        y = A*sin(w*t+phi);     % excitation this input creates a current in the voice coil
        
        cs=[cos(w*t) sin(w*t)];	% build matrix in one command to find the Sine/Cose coefficients
        %fit sine
        a = cs\y;               % solve least squares problem exploiting Matlab's \
        yfit = cs*a;            % predict y using fitted model
        
        figure(1)
        h = plot(t,y,'o',t,yfit,t,yfit-y);
        set(h,'LineWidth',2);
        hl = legend('$y_{measured}$','$y_{fitted}$','Error');
        set(hl,'interpreter','latex','FontSize',18)
        xlabel('time[sec]','interpreter','latex','FontSize',18)
        set(gcf,'Color','w')
        set(gca,'FontSize',14)
        xlim([t(1) t(end)])
        drawnow, shg
        shg
        
        A_fitted = norm(a); phi_fitted = atan2(a(1),a(2))*180/pi;
        if phi_fitted<0, phi_fitted = phi_fitted+360; end % just for display
        
        title(sprintf('$A_{fit}=%.4f$,  $\\phi_{fit}=%.1f$',A_fitted ,phi_fitted),'interpreter','latex','FontSize',18)
        
    case 'Stroboscope demo'
        
        x=0:.01:5;FS=500;
        t=0:1/FS:50;f=60;T=1/f;Nsam=T*FS;NsamNoStob=1/25/0.01;
        y=sin(2*pi*f*t)'*sin(2*pi*x);
        for k=1:length(t)/5
            subplot(3,1,1);plot(x,y(fix(k*NsamNoStob),:));ylim([-2 2]);title('No stroboscope');
            subplot(3,1,2);plot(x,y(fix(k*Nsam/1.03),:));ylim([-2 2]);title('Stroboscope freq is 1.03 times bigger');
            subplot(3,1,3);plot(x,y(fix(k*Nsam/1.07),:));ylim([-2 2]);title('Stroboscope freq is 1.07 times bigger');
            drawnow;pause(.1);
        end
        
end