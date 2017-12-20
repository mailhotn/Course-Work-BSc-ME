%% Q2 Z-N
P = zpk([],-1,1);
rlocus(P)

%% Q3 Z-N
P = zpk([],[-1 -1 -1],1);
k = linspace(0,100,1001);
rlocus(P,k)
ku = 8;
Gu = feedback(ku*P,1);
[mag,~,w] = bode(Gu);
[~,ind] = max(mag);
Tu = 2*pi/w(ind);
C1 = pidstd(0.45*ku,Tu/1.2);
C2 = pidstd(0.6*ku,Tu/2,Tu/8);

%% Q5 Sim
k_p = C1.Kp;
k_i = 1/C1.Ti;
s = tf('s');
T_ABB = P*k_i*k_p/(s*(P*k_i/s*k_p+k_p*P+1));

step(feedback(C1*P,1),T_ABB)
legend('Classic PI','ABB PI')