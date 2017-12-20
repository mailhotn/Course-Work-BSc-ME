%% Model
clc;h=.6e-3;b=30e-3;L=128.5e-3;R=40e-3;m=300e-3;E=200e9;I=b*h^3/12;
M=m*eye(4);k1=12*E*I/L^3;k2=E*I/R^3/(pi/2-4/pi);
K = [2*k1+k2    -k2           0            0
     -k2      2*k1+2*k2      -k2          0
     0          -k2        2*k1+2*k2      -k2
     0           0           -k2        2*k1+k2];
C = 0.1*eye(4);

A = [zeros(4) eye(4); -M\[K C]];
B = [zeros(4);inv(M)];
c_sys =[eye(4) zeros(4)];
D = zeros(4);
sys = ss(A,B,c_sys,D); 

%% Chirp
Fs  = 500; %sampling frequency
f0  = 1;
f1  = 60;
t1  = 15;
t   = 0:1/Fs:t1; %time base - upto 1 second
chu = chirp(t,f0,t1,f1);

%% Simulation
u = [chu; zeros(1,length(t)); zeros(1,length(t)); zeros(1,length(t))];
y = lsim(sys,u,t);
figure(1)
plot(t,y)
xlabel('t [sec]')
ylabel('Displacement [mm]')
legend('m_1','m_2','m_3','m_4')