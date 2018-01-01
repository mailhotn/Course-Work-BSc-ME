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
for ii = 1:epochs
    ni = X_train*w;
    e = T_train - ni;
    dw = eta*X_train.'*e;
    Jp = [Jp, 0.5*e.'*e]; %#ok
    w = w + dw;
%     dJ = Jp(end-1)-Jp(end);
end

figure(1)
plot(1:length(Jp),Jp)
ylabel('Mean Square Error')
xlabel('Epoch')

ni = X_test*w;
u = sign(ni);
e = T_test-u;
N_errors_test = e.'*e/4

ni = [ones(347,1) X]*w;
u = sign(ni);
e = T-u;
N_errors_all = e.'*e/4

%% Q2 Reduced Order
load('home21_dat.mat');

data = data(randperm(size(data,1)),:);
T = data(:,1);
T(find(T==0)) = -1;
X = data(:,2:11);
X(:,[1, 2, 5, 6, 7, 8, 9]) = []; % Only take 3, 4, 10
T_train = T(1:300);
X_train = [ones(300,1),X(1:300,:)];
T_test = T(301:347);
X_test = [ones(47,1),X(301:347,:)];

w  = 2*rand(size(X_train,2),1)-1;
Jp = Inf;
dJ = Inf;
eta = 1e-6;
epochs = 10000;
for ii = 1:epochs
    ni = X_train*w;
    e = T_train - ni;
    dw = eta*X_train.'*e;
    Jp = [Jp, 0.5*e.'*e]; %#ok
    w = w + dw;
    dJ = Jp(end-1)-Jp(end);
end

figure(1)
plot(1:length(Jp),Jp)
ylabel('Mean Square Error')
xlabel('Epoch')
ni = X_test*w;
u = sign(ni);
e = T_test-u;
N_errors_test = e.'*e/4

ni = [ones(347,1) X]*w;
u = sign(ni);
e = T-u;
N_errors_all = e.'*e/4

%% Reduced Order Scatter
load('home21_dat.mat');
T = data(:,1);
T(find(T==0)) = -1;
X = data(:,2:11);
X(:,[1, 2, 5, 6, 7, 8, 9]) = [];
plot3(X(find(T==1),1),X(find(T==1),2),X(find(T==1),3),'bx',X(find(T==-1),1),X(find(T==-1),2),X(find(T==-1),3),'r+')
xlabel('Feature 3')
ylabel('Feature 4')
zlabel('Feature 10')