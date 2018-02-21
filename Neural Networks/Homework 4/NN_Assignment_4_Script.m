%% Question 1 - Classification in half-rings
[A_train, B_train] = RandInRing(10,6,-4,1000);
X_train = [A_train; B_train];
T_train = [ones(1,1000), zeros(1,1000);
           zeros(1,1000), ones(1,1000)];     
[A_val, B_val] = RandInRing(10,6,-4,500);
X_val = [A_val; B_val];
T_val = [ones(1,500), zeros(1,500);
         zeros(1,500),ones(1,500)];
[A_test, B_test] = RandInRing(10,6,-4,1000);




X_test = [A_test; B_test];
T_test = [ones(1,1000), zeros(1,1000);
          zeros(1,1000), ones(1,1000)];
X = [X_train;
     X_val;
     X_test].';
T = [T_train, T_val, T_test];

sse_tr = 0*(2:2:20);
sse_val = 0*(2:2:20);
for ii = 1:10
    net = patternnet(2*ii);
    net.trainFcn = 'trainlm';
    net.performFcn = 'sse';
    net.divideFcn = 'divideind';
    net.divideParam.trainInd = 1:2000;
    net.divideParam.valInd   = 2001:3000;
    net.divideParam.testInd  = 3001:5000;
    [net,tr] = train(net,X,T); %#ok
    sse_tr(ii) = tr.perf(end);
    sse_val(ii) = tr.vperf(end);
end

figure(1)
plot(2:2:20,sse_tr,'o',2:2:20,sse_val,'o')
legend('Training','Validation')
ylabel('Performance (SSE)')
xlabel('Hidden Layer Neurons')

%% 4 is enough
% Train
net = patternnet(3);
net.trainFcn = 'trainlm';
net.performFcn = 'sse';
net.divideFcn = 'divideind';
net.divideParam.trainInd = 1:2000;
net.divideParam.valInd   = 2001:3000;
net.divideParam.testInd  = 3001:5000;
[net,tr] = train(net,X,T);
% Decision Boundary
x = -20:1:30;
y = -20:1:20;
[xx, yy] = meshgrid(x,y);
zz = 0*xx;
for ii = 1:length(y)
    for jj = 1:length(x)
        activation = sim(net,[xx(ii,jj) yy(ii,jj)].');
        if activation(1) > activation(2)
            zz(ii,jj) = 1;
        else
            zz(ii,jj) = 0;
        end
    end
end
% Plot stuff
figure(2)
axis([-20 30 -20 20])
axis equal
grid on
colormap([0 0 0; 0 0 0])
contour(xx,yy,zz,1,'Linewidth',2)
hold on
plot(A_train(:,1),A_train(:,2),'bo',B_train(:,1),B_train(:,2),'rx')
hold off

%% Question 2 - PCA
load('home21_dat.mat');
samples = data(:,2:end);
T       = [data(:,1).';
    1 - data(:,1).'];
A       = zscore(samples);

% PCA using SVD
[U,S,~] = svd(A,'econ');
s       = diag(S);

for ii = 1:10
    sums(ii) = sum(s(1:ii).^2); %#ok
end
figure(); plot(1:10, s.^2/(sum(s.^2)),'--o',1:10, sums/sum(s.^2),'--o');
xlabel('\sigma_{i}'); ylabel('%Variance and accumulated %Variance');
legend('%Variance','accumulated %Variance');



% classification using some of the data
[coeff, score, latent] = pca(A);
for ii = 1:10
    net      = patternnet(2);
    net.performFcn = 'sse';
    [net,tr] = train(net,(score(:,1:ii)).',T);
    best(ii) = tr.best_vperf; %#ok
end
figure; plot(1:10, best, 'o');
xlabel('PC index'); ylabel('Validation Group best SSE');