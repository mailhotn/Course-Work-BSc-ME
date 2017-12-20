% Preparation report1 Shai 16-Oct-2016
% clc;h=.5e-3;b=30e-3;L=128.5e-3;R=40e-3;m=200e-3;E=200e9;I=b*h^3/12;
clc;h=.6e-3;b=30e-3;L=128.5e-3;R=40e-3;m=300e-3;E=200e9;I=b*h^3/12;
M = m*eye(4);k1=12*E*I/L^3;k2=E*I/R^3/(pi/2-4/pi);
K = [2*k1+k2    -k2           0            0
     -k2      2*k1+2*k2      -k2          0
     0          -k2        2*k1+2*k2      -k2
     0           0           -k2        2*k1+k2];
%q1_b
[V, W2]=eig(K,M);wn = sqrt(diag(W2))/2/pi;
%q1_c
C = 0.1*eye(4);
% State-Space
A = [zeros(4) eye(4); -M\[K C]];
B = [zeros(4);inv(M)];
c_sys =[eye(4) zeros(4)];
D = zeros(4);
sys = ss(A,B,c_sys,D); 
% set(cstprefs.tbxprefs,'FrequencyUnits','Hz','FreqScale','linear');
opts = bodeoptions('cstprefs');
opts.FreqUnits = 'Hz';
opts.FreqScale = 'linear';
figure(1);bodeplot(sys(1,1),sys(2,1),sys(3,1),sys(4,1),opts);xlim([5 45]);
legend('m_1','m_2','m_3','m_4')
figure(2);bode(sys(2,1));xlim([5 100]);
figure(3);bode(sys(3,1));xlim([5 100]);
figure(4);bode(sys(4,1));xlim([5 100]);
%q2_a
figure(1);clf;step(sys(:,1));[y,t]=step(sys);
%q2_b
Fs = 1/(t(2)-t(1));
[X_1 f]=FourierSeriesDFT(y(:,1,1),Fs);
figure(1);clf;semilogy(f, abs(X_1));xlim([0 60]);
xlabel 'frequency Hz'
ylabel |X_1(\omega)| 
[X_2 f]=FourierSeriesDFT(y(:,1,2),Fs);
figure(2);clf;
semilogy(f, abs(X_2));xlim([0 70]);
xlabel 'frequency Hz'
ylabel |X_2(\omega)| 
[X_3 f]=FourierSeriesDFT(y(:,1,3),Fs);
figure(3);clf;
semilogy(f, abs(X_3));xlim([0 70]);
xlabel 'frequency Hz'
ylabel |X_3(\omega)| 
[X_4 f]=FourierSeriesDFT(y(:,1,4),Fs);
figure(4);clf;
semilogy(f, abs(X_4));xlim([0 70]);
xlabel 'frequency Hz'
ylabel |X_4(\omega)|
%q3_b
t_input_end=30;
t_no_in_end=t_input_end+20;
dt=t_no_in_end/15000;
t_sim=0:dt:t_no_in_end;
t_input=0:dt:t_input_end;
t_no_in=(t_input_end+dt):dt:t_no_in_end;

for ii = 1:4
    close all
    % 1st freq
    u_w=[sin(wn(ii)*2*pi*t_input') zeros(size(t_input',1),3);
        zeros(size(t_no_in',1),1) zeros(size(t_no_in',1),3)];
    figure(1);clf;
    lsim(sys,u_w,t_sim')
    y_lsim_w1=lsim(sys,u_w,t_sim');
    damping_coeff(t_sim' ,y_lsim_w1, wn(ii))
end
 
%q4
figure(1);clf;
t_q5=0:30/10000:30;
in1=[sin(3*t_q5') zeros(size(t_q5',1),3)];
in2=[sin(7*t_q5') zeros(size(t_q5',1),3)];
figure(2);clf;
lsim(sys,in1,t_q5')
y1=lsim(sys,in1,t_q5');
figure(3);clf;
lsim(sys,in2,t_q5')
y2=lsim(sys,in2,t_q5');
figure(4);clf;
lsim(sys,(2*in1+8*in2),t_q5')
ysum=lsim(sys,(2*in1+8*in2),t_q5');
diff=ysum-(2*y1+8*y2);
figure(5);clf;
plot(t_q5, diff)
xlabel 'time sec'
ylabel 'difference'