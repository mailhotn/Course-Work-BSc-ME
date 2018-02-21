%% Question 1
%Part a
P = tf(1,[0.064, 0.02, 1.176]);
Ci = tf([1, 4.95, 2.2],[0.025,1,0]);
% sisotool(Ci*P)
C = 0.01586*Ci;
figure(1)
bode(C*P)
figure(2)
step(feedback(P,C))

%Part b
figure(1)
bode(feedback(C*P,1))
fb = 4.42/2*pi;
ts = 1/(2.2*fb);

%Part c
% Equivalent Delay
h_eq = ts/2;
P1 = P;
P1.inputdelay = h_eq;
figure(1)
step(feedback(P,C),'-b',feedback(P1,C),'-r',0:ts:180)

% Discrete Time
P2 = c2d(P,ts);
Cd = c2d(C,ts,'tustin');
figure(2)
step(feedback(P,C),'-b',feedback(P2,Cd),'-r',0:ts:180)

% Multi-Rate Simulation
load_system('Q1_Part_c');
set_param('Q1_Part_c','StopTime','180','AbsTol','1e-10','RelTol','1e-10');
sim('Q1_Part_c');
figure(3)
step(feedback(P,C),'-b',0:ts:180)
hold on
plot(Output.time, Output.signals.values,'-r');
hold off

%Part d Anti - Aliasing
wn = pi/ts;
F = zpk(z,p,k);
[z,p,k] = cheby2(2,20,wn,'s');
figure(1)
step(feedback(P,C),feedback(P,C*F),0:ts/10:180)
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
figure(); bode(minreal(L));
Gyr = feedback(C2*P2*feedback(P1*C1,1),1);
step(Gyr);
load_system('Q2_a1');
set_param('Q2_a1','StopTime','1','AbsTol','1e-10','RelTol','1e-10');
sim('Q2_a1');
figure(); plot(cscdOut.time, cscdOut.signals.values);

% Cascade design simulation
load_system('Q2_b1');
set_param('Q2_b1','StopTime','1','AbsTol','1e-10','RelTol','1e-10');
sim('Q2_b1');
figure(); plot(cscdOut.time, cscdOut.signals(1).values,...
    cscdOut.time, cscdOut.signals(2).values);

% Using position measurement to estimate angular velocity
e = 1e-1;
C = Cpi*Cp;
H = 1 + 1/Cp*tf([1 0],[e 1]);
H = minreal(H);
load_system('Q2_b2');
set_param('Q2_b2','StopTime','2','AbsTol','1e-10','RelTol','1e-10');
sim('Q2_b2');
figure(); plot(cscdOut.time, cscdOut.signals(1).values,...
    cscdOut.time, cscdOut.signals(2).values);

% Controller designed directly
load_system('Q2_Part_b3');
load('Q2_a3_C.mat');
sim('Q2_Part_b3');
set_param('Q2_Part_b3','StopTime','1','AbsTol','1e-8','RelTol','1e-8');
figure()
plot(Sim_Out_y.time, Sim_Out_y.signals.values,...
    Sim_Out_u.time, Sim_Out_u.signals.values);