
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

load('C.mat');
P = tf(1/(N*Kt), [R*(Jl + Jm)/(N^2*Kt^2), R*c/(N^2*Kt^2) + 1]);
L = C*P;
figure();hold on
step(feedback(L,1)*pi/2,10);
for ii = 1:2  
    A = [0 1 0 0;
        -K(ii)/Jm -Kt^2/(R*Jm) N*K(ii)/Jm 0;
        0 0 0 1;
        N*K(ii)/Jl 0 -N^2*K(ii)/Jl -c/Jl];
    B = [0 Kt/(R*Jm) 0 0].';
    C = [0 0 0 1]; % here C gives angle
    D = 0;
    P = ss(A,B,C,D);
    load('C.mat');
    step(pi/2*feedback(C*P,1),10);
end
load_system('DCMotorLoopwBacklashAngularVel');
set_param('DCMotorLoopwBacklashAngularVel','StopTime','10','AbsTol','1e-10','RelTol','1e-10');
t = {'r', 'b'};
for ii = 1:2  
    load('C.mat');
    sim('DCMotorLoopwBacklashAngularVel');
    plot(Output.time, Output.signals.values,'color', t{ii});   
end
hold off
xlabel('Time [sec]'); ylabel('d\theta_{l} [rad/sec]');




load('C_PI.mat');
P = tf(1/(N*Kt), [R*(Jl + Jm)/(N^2*Kt^2), R*c/(N^2*Kt^2) + 1, 0]);
L = C*P;
figure();
hold on
step(feedback(L,1)*pi/6,60);
for ii = 1:2  
    A = [0 1 0 0;
        -K(ii)/Jm -Kt^2/(R*Jm) N*K(ii)/Jm 0;
        0 0 0 1;
        N*K(ii)/Jl 0 -N^2*K(ii)/Jl -c/Jl];
    B = [0 Kt/(R*Jm) 0 0].';
    C = [0 0 1 0]; % here C gives angle
    D = 0;
    P = ss(A,B,C,D);
    load('C_PI.mat');
    step(pi/6*feedback(C*P,1),60);
end

load_system('DCMotorLoopwBacklashAngle');
set_param('DCMotorLoopwBacklashAngle','StopTime','60','AbsTol','1e-10','RelTol','1e-10');
t = {'r', 'b'};
for ii = 1:2  
    load('C_PI.mat');
    sim('DCMotorLoopwBacklashAngle');
    plot(Output.time, Output.signals.values,'color',t{ii});   
end
xlabel('Time [sec]'); ylabel('\theta_{l} [rad]');
hold off
