%% Q1 - Fly
m = 0.001;
k = 0.2e6;
c = 2;
a = 0.1e-6;
V_max = 100;
L = 100e-6; %#ok
l = 280e-6;
th_max = 10/180*pi;
c_a = 9e-5;

% Calculate needed layers
n = ceil(l*th_max/a/V_max);

% Build Model
A = [0 m zeros(1,8);
    -2*k -2*c k c zeros(1,6);
     zeros(1,3) m zeros(1,6);
     k c -2*k -2*c k c zeros(1,4);
     zeros(1,5) m zeros(1,4);
     zeros(1,2) k c -2*k -2*c k c zeros(1,2);
     zeros(1,7) m zeros(1,2);
     zeros(1,4) k c -2*k -2*c k c;
     zeros(1,9) m;
     zeros(1,6) k c -k -c]/m;
B = [zeros(9,2);
     k*a/m 1/m]; % u_1 = V, u_2 = F_ext
C = [zeros(1,8) 1 0]; % Sum of displacements
D = 0;
sys = ss(A,B(:,1),C,D);
figure(1)
bode(sys)

% Angle Out
C2 = -C./l;
B2 = [zeros(9,1);
     k*a/m];
A_add = -2*c_a/l^2*[zeros(1,9) 1];
A2 = A + B(:,2)*A_add;
D2 = 0;
sys2 = ss(A2,B2,C2,D2);
figure(2)
bode(sys2)

% Required Amplitude
mag = bode(sys2,600*2*pi);
Amp = th_max/mag;
%% Q2 - Camera on Cart
dx_max = 30*1000/3600;
ddx_max = 25;
J = 0.07;
m = 0.8;
l = 0.3;
c_a = 0.1;
R = 0.76;
K_t = 0.46;
V_max = 20;
D_v = 1e-2;
T_dist = c_a*dx_max^2 + m*l*ddx_max;
T_stall = K_t*V_max/R;
w_free = V_max/K_t;
% plot Graphs
X = [zeros(4,1) w_free*[1/1 1/3 1/5 1/7].'].';
Y = [T_stall*[1*0.98;3*0.95;5*0.9;7*0.88] zeros(4,1)].';
plot(X,Y,0,T_dist,'b*')

%% Q3 - Servo motor control

Kt   = 0.03;       % [Nm/A]
R    = 0.1;        % [Ohm]
L    = 0; %#ok     % [H]
Vmax = 1;          % [V]
Jm   = 0.001;      % [Kgm^2] 
Jl   = 0.02;       % [Kgm^2]
c    = 0.001;      % [Nmsec/rad]
N    = 3;
K    = [0.1 30];   % [Nm/rad]
d    = 5*pi/180;   % [rad]

%-------------------------------------------------------------------------%
%                           rigid transmission
%-------------------------------------------------------------------------%

% angular velocity control
P = tf(1/(N*Kt), [R*(Jl + Jm)/(N^2*Kt^2), R*c/(N^2*Kt^2) + 1]);
load('C.mat');
L = C*P;
% sisotool(C*P); % design controller (OS < 20 %)
[Gm, Pm, Wgm, Wpm] = margin(L) %#ok
S            = stepinfo(feedback(L,1)) %#ok
S_controlSig = stepinfo(feedback(C,P)); 
maxReference = 1/(S_controlSig.Peak)*180/pi %#ok %[deg/sec] 
figure(); step(feedback(L,1)*pi/2);

% angle control
P = tf(1/(N*Kt), [R*(Jl + Jm)/(N^2*Kt^2), R*c/(N^2*Kt^2) + 1, 0]);
% sisotool(P); % design controller (PM > 70)
load('C_PI.mat');
L = C*P;
[Gm, Pm, Wgm, Wpm] = margin(L) %#ok
S            = stepinfo(feedback(L,1)) %#ok
S_controlSig = stepinfo(feedback(C,P)); 
maxReference = 1/(S_controlSig.Peak)*180/pi %#ok %[deg/sec] 
figure(); step(feedback(L,1)*pi/6);

%-------------------------------------------------------------------------%
%                           flexible transmission
%-------------------------------------------------------------------------%

% State-space model of the plant (see report for derivation)
% X = [th_m dth_m th_l dth_l]^T
for ii = 1:2  
    A = [0 1 0 0;
        -K(ii)/Jm -Kt^2/(R*Jm) N*K(ii)/Jm 0;
        0 0 0 1;
        N*K(ii)/Jl 0 -N^2*K(ii)/Jl -c/Jl];
    B = [0 Kt/(R*Jm) 0 0].';
    C = [0 0 0 1]; % here C gives angular velocity
    D = 0;
    P = ss(A,B,C,D);
    load('C.mat');
    [Gm, Pm, Wgm, Wpm] = margin(C*P) %#ok
    figure(); step(pi/2*feedback(C*P,1));
    A = [0 1 0 0;
        -K(ii)/Jm -Kt^2/(R*Jm) N*K(ii)/Jm 0;
        0 0 0 1;
        N*K(ii)/Jl 0 -N^2*K(ii)/Jl -c/Jl];
    B = [0 Kt/(R*Jm) 0 0].';
    C = [0 0 1 0]; % here C gives angle
    D = 0;
    P = ss(A,B,C,D);
    load('C_PI.mat');
    [Gm, Pm, Wgm, Wpm] = margin(C*P) %#ok
    figure(); step(pi/6*feedback(C*P,1));
end

load_system('DCMotorLoopwBacklashAngle');
set_param('DCMotorLoopwBacklashAngle','StopTime','60','AbsTol','1e-10','RelTol','1e-10');
for ii = 1:2  
    load('C_PI.mat');
    sim('DCMotorLoopwBacklashAngle');
    figure(); plot(Output.time, Output.signals.values);
    xlabel('Time [sec]'); ylabel('\theta_{l} [rad]');
end

load_system('DCMotorLoopwBacklashAngularVel');
set_param('DCMotorLoopwBacklashAngularVel','StopTime','10','AbsTol','1e-10','RelTol','1e-10');
for ii = 1:2  
    load('C.mat');
    sim('DCMotorLoopwBacklashAngularVel');
    figure(); plot(Output.time, Output.signals.values);
    xlabel('Time [sec]'); ylabel('d\theta_{l} [rad/sec]');
end