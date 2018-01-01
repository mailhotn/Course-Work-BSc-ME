%% Q1 Part 1
load('home11_dat.mat');

N = 3;
n = con2seq(n);
m = con2seq(m);
net = newlin([-30, 30], 1, 0:(N-1));
net.IW{1,1} = rand(1,N);
net.b{1}    = rand(1);
[net,Y,e,Pf] = adapt(net,n,m);
sound(cat(2,e{:}));
figure(1)
plot(cat(2,e{:}));

%% Q1 Part 2
load('home12_dat.mat');

N = 4;
n = con2seq(n);
m = con2seq(m);
net = newlin([-30, 30], 1, 0:(N-1));
net.IW{1,1} = rand(1,N);
net.b{1}    = rand(1);
[net,Y,e,Pf] = adapt(net,n,m);
sound(cat(2,e{:}));
figure(2)
plot(cat(2,e{:}));

%% Q2
load('home21_dat.mat');

data = data(randperm(size(data,1)),:);
T = data(:,1);
T(find(T==0)) = -1;
X = data(:,2:11);
T_train = T(1:300);
X_train = [ones(300,1),X(1:300,:)];
T_test = T(301:347);
X_test = [ones(47,1),X(301:347,:)];

w  = 2*rand(11,1)-1;
Jp = Inf;
dJ = Inf;
eta = 1e-9;
epochs = 100000;
while abs(dJ) > 1e-4
    ni = X_train*w;
    e = T_train - ni;
    dw = eta*X_train.'*e;
    Jp = [Jp, 0.5*e.'*e]; %#ok
    w = w + dw;
    dJ = Jp(end-1)-Jp(end);
end

figure(1)
plot(1:length(Jp),Jp)
ni = X_test*w;
u = sign(ni);
e = T_test-u;
J_test = e.'*e/4

%% Reduced Order
load('home21_dat.mat');

data = data(randperm(size(data,1)),:);
T = data(:,1);
T(find(T==0)) = -1;
X = data(:,2:11);
X(:,[2, 4, 5, 6, 7, 8, 9]) = [];
T_train = T(1:300);
X_train = [ones(300,1),X(1:300,:)];
T_test = T(301:347);
X_test = [ones(47,1),X(301:347,:)];

w  = 2*rand(11,1)-1;
Jp = Inf;
dJ = Inf;
eta = 1e-9;
epochs = 100000;
while abs(dJ) > 1e-4
    ni = X_train*w;
    e = T_train - ni;
    dw = eta*X_train.'*e;
    Jp = [Jp, 0.5*e.'*e]; %#ok
    w = w + dw;
    dJ = Jp(end-1)-Jp(end);
end

figure(1)
plot(1:length(Jp),Jp)
ni = X_test*w;
u = sign(ni);
e = T_test-u;
J_test = e.'*e/4