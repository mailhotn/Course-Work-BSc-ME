%% Question 1
%Part a
P = tf(1,[0.064, 0.02, -1.176]);
Ci = tf([1, 4.95, 2.2],[0.025,1,0]);
% sisotool(Ci*P)
C = Ci;
figure(1)
bode(C*P)
figure(2)
step(feedback(P,C))

%Part b
figure(1)
bode(feedback(C*P,1))
fb = 24.6/2*pi;
ts = 1/(2.2*fb);

%Part c
% Equivalent Delay
h_eq = ts/2;
P1 = P;
P1.inputdelay = h_eq;
figure(1)
step(feedback(P,C),'-b',feedback(P1,C),'-r');

% Discrete Time
P2 = c2d(P,ts);
Cd = c2d(C,ts,'tustin');
figure(2)
step(feedback(P,C),'-b',feedback(P2,Cd),'-r',0:ts:14);

candidate_ts = [];
for ts_s = linspace(9*ts, 12*ts, 1000)
    P2_s = c2d(P,ts_s);
    Cd_s = c2d(C,ts_s,'tustin');
    info = stepinfo(feedback(P2_s,Cd_s));
    if info.Peak > 10
        disp('Found one');
        candidate_ts = [candidate_ts, ts_s]; %#ok
        step(feedback(P2_s,Cd_s));
        pause;
    end
end


% Multi-Rate Simulation
load_system('Q1_Part_c');
set_param('Q1_Part_c','StopTime','14','AbsTol','1e-10','RelTol','1e-10');
sim('Q1_Part_c');
figure(3)
step(feedback(P,C),'-b',0:ts:14)
hold on
plot(Output.time, Output.signals.values,'-r');
hold off

candidate_ts = [];
for ts_s = linspace(10.5*ts, 11*ts, 500)
    Cd_s = c2d(C,ts_s,'tustin');
    sim('Q1_Part_c3');
    info = stepinfo(Output.signals.values, Output.time);
    if info.Peak > 10
        disp('Found one');
        candidate_ts = [candidate_ts, ts_s]; %#ok
        plot(Output.time, Output.signals.values,'-r');
        pause;
    end
end

%Part d Anti - Aliasing
wn = pi/ts;
[z,p,k] = cheby2(2,20,wn,'s');
F = zpk(z,p,k);
load_system('Q1_Part_d');
set_param('Q1_Part_d','StopTime','14','AbsTol','1e-10','RelTol','1e-10');
sim('Q1_Part_d');
figure();
step(feedback(P,C),'-b',0:ts/100:14)
hold on
plot(Output.time, Output.signals.values,'-r');
hold off

%% Question 2
enc_res = 0.005; % encoder resolution
Kt = 0.46;  % [Nm/A]
R  = 0.76;  % [Ohm]
Jl = 0.025; % [kgm^2]
Ke = Kt;
tau_m = Jl*R/(Kt*Ke);
P  = 1/Ke*tf(1, [tau_m 1]);
load('Q2Cp.mat');

wb  = 150; % [rad/sec]
Cpi = P^-1*wb*tf(1, [1 0]);
C1 = Cpi;
C2 = Cp;
P1 = P; 
P2 = tf(1,[1 0]);
L  = C1*C2*P2*P1 + C1*P1;
[Gm_s, Pm_s, Wgm_s, Wpm_s] = margin(minreal(L));
% Frequency
figure(); bode(minreal(L));
title('')
% Reference Step
Gyr = feedback(C2*P2*feedback(P1*C1,1),1);
figure(); step(Gyr);
title('')
% Disturbance Step
load_system('Q2_a1');
set_param('Q2_a1','StopTime','1','AbsTol','1e-10','RelTol','1e-10');
sim('Q2_a1');
figure(); plot(cscdOut.time, cscdOut.signals.values);
xlabel('Time [sec]')
ylabel('Angle [rad]')

% Cascade design simulation
load_system('Q2_b1');
set_param('Q2_b1','StopTime','0.4','AbsTol','1e-10','RelTol','1e-10');
sim('Q2_b1');
figure(); plot(cscdOut.time, cscdOut.signals(1).values,...
    cscdOut.time, cscdOut.signals(2).values);
xlabel('Time [sec]')
ylabel('Amplitude')
legend('Response','Control Effort')

% Using position measurement to estimate angular velocity
e = 1e-3;
C = Cpi*Cp;
H = 1 + 1/Cp*tf([1 0],[e 1]);
H = minreal(H);
% Frequency
figure(); bode(minreal(C*P1*P2*H));
title('')
% Step
figure(); step(feedback(C*P1*P2,H));
title('')
% Disturbance Step
figure(); step(feedback(P1*P2,C*H));
title('')
% Simulation
load_system('Q2_b2');
set_param('Q2_b2','StopTime','2','AbsTol','1e-7','RelTol','1e-7');
sim('Q2_b2');
figure(); plot(cscdOut.time, cscdOut.signals(1).values,...
    cscdOut.time, cscdOut.signals(2).values);
xlabel('Time [sec]')
ylabel('Amplitude')
legend('Response','Control Effort')

% Controller designed directly
load('Q2_a3_C.mat');
% Frequency
figure(); bode(minreal(C*P1*P2));
title('')
% Step
figure(); step(feedback(C*P1*P2,1));
title('')
% Disturbance Step
figure(); step(feedback(P1*P2,C));
title('')
% Simulation
load_system('Q2_Part_b3');
sim('Q2_Part_b3');
set_param('Q2_Part_b3','StopTime','1','AbsTol','1e-8','RelTol','1e-8');
figure()
plot(Sim_Out_y.time, Sim_Out_y.signals.values,...
    Sim_Out_u.time, Sim_Out_u.signals.values);
xlabel('Time [sec]')
ylabel('Amplitude')
legend('Response','Control Effort')

%% Question 3
% Given
K_Fa = 2.16e5;
K_Fd = -5.4e4;
K_Ta = -2.7e5;
K_Tq = -2.25e3;
K_Td = 3.24e5;
I    = 800;
m    = 1e3;
V    = 300;
wn   = 500;
z    = 1.2;
H = tf(wn^2,[1, 2*z*wn, wn^2]);

% Part a - Controller Design
A = [-K_Fa/(m*V) 1;
    K_Ta/I K_Tq/I];
B = [-K_Fd/(m*V) K_Td/I].';
C = [K_Fa/m 0;
    0   1];
D = [K_Fd/m 0].';

P = ss(A,B,C,D);
Pa = tf(P(1));
Pq = tf(P(2));
% sisotool
% K1 = C2;
% K2 =  C1.K/K1;
load('Q3_a_K123.mat')
Gd = minreal(K1*K2/(tf('s') + K1*Pq*H*(K2+tf('s'))));
% sisotool(Pa*Gd)
L2r = Pa*Gd*K3;
L = K1*K2*K3/tf('s')*Pa + K1*(1+K2/tf('s'))*Pq;
figure()
bode(L2r,L)
legend('L2r','L')

load_system('Q3_a_15');
set_param('Q3_a_15','StopTime','1','AbsTol','1e-8','RelTol','1e-8');
sim('Q3_a_15');
figure()
plot(Sim_Out.time, Sim_Out.signals.values);
hold on
info = stepinfo(Sim_Out.signals.values,Sim_Out.time) %#ok

% Part b - Moved Sensor
L = -0.7;
A = [-K_Fa/(m*V) 1;
    K_Ta/I K_Tq/I];
B = [-K_Fd/(m*V) K_Td/I].';
C = [K_Fa/m+L*K_Ta/I L*K_Tq/I;
    0   1];
D = [K_Fd/m+L*K_Td/I 0].';

P = ss(A,B,C,D);
Pa_m = tf(P(1));
Pq = tf(P(2));

% Controller Correction
K1_b = K1*(1-K2*K3*L);
K2_b = K2/(1-K2*K3*L);

% Simulation with Controller correction
load_system('Q3_b_15');
set_param('Q3_b_15','StopTime','1','AbsTol','1e-8','RelTol','1e-8');
sim('Q3_b_15');
plot(Sim_Out.time, Sim_Out.signals.values);
xlabel('Time [sec]')
ylabel('Acceleration [m/s^2]')
legend('Sensor at COM','Sensor not at COM')
hold off
info = stepinfo(Sim_Out.signals.values,Sim_Out.time) %#ok

% Stability margins with controller correction
Gd = minreal(K1_b*K2_b/(tf('s') + K1_b*Pq*H*(K2_b+tf('s'))));
L2r = Pa_m*Gd*K3;
L = K1_b*K2_b*K3/tf('s')*Pa_m + K1_b*(1+K2_b/tf('s'))*Pq;
figure()
bode(L2r,L)
legend('L2r','L')

% System Response Without Correction
load_system('Q3_b1_15');
set_param('Q3_b1_15','StopTime','1','AbsTol','1e-8','RelTol','1e-8');
sim('Q3_b1_15');
figure()
plot(Sim_Out.time, Sim_Out.signals.values);
xlabel('Time [sec]')
ylabel('Acceleration [m/s^2]')
info = stepinfo(Sim_Out.signals.values,Sim_Out.time) %#ok