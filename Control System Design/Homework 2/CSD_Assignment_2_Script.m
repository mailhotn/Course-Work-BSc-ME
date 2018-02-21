%% Question 1 - Part a
P = tf(1,[0.064, 0.02, 1.176]);
Ci = tf([1, 4.95, 2.2],[0.025,1,0]);
% sisotool(Ci*P)
C = 0.01586*Ci;
figure(1)
bode(C*P)
figure(2)
step(feedback(P,C))

%% Question 1 - Part b
figure(1)
bode(feedback(C*P,1))
fb = 25.8/2*pi;
ts = 1/(2.2*fb);

%% Question 1 - Part c
% Equivalent Delay
h_eq = ts/2;
P1 = P;
P1.inputdelay = h_eq;
figure(1)
step(feedback(P,C),'-b',feedback(P1,C),'-r')

% Discrete Time
P2 = c2d(P,ts);
Cd = c2d(C,ts,'tustin');
figure(2)
step(feedback(P,C),'-b',feedback(P2,Cd),'-r')

% Multi-Rate Simulation
load_system('Q2_Part_c');
set_param('Q2_Part_c','StopTime','180','AbsTol','1e-10','RelTol','1e-10');
sim('Q2_Part_c');
figure(3)
hold on
step(feedback(P,C),'-b')
plot(Output.time, Output.signals.values,'-r');
hold off