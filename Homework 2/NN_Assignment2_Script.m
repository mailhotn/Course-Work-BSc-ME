%% Q2

% init parameters
mu = [0.5 0.5];
Sigma = [0.1 0.05;
    0.05 0.1];
x1 = -1:.03:2; x2 = x1;
[X1, X2] = meshgrid(x1,x2);

% calculate the pdf of each point, then plot
Fw2 = mvnpdf([X1(:), X2(:)], mu, Sigma);
Fw2 = reshape(Fw2, length(x2), length(x1));
figure()
contour(x1,x2,Fw2,[.0001 .001 .01 1]);
line([0 0 1 1 0], [0 1 1 0 0],'linestyle','--','color','k');
text(0.5, 0.5, '\omega_{2}'); text(0.75, 0.25,'\omega_{1}');
text(1.2, 0.5, '\omega_{2}');
xlabel('x_{1}'); ylabel('x_{2}');
axis equal
hold off
figure()
surf(x1,x2,Fw2);
hold on
v = [0 0 1; 1 0 1; 1 1 1; 0 1 1];
f = 1:4;
patch('Faces',f,'Vertices',v,'EdgeColor','b',...
    'FaceColor','none','linewidth',2);
hold off
xlabel('x_{1}'); ylabel('x_{2}'); zlabel('Probability Density');

% Whitening transformation
[V, D] = eig(Sigma);
Aw = V*D^(-1/2);
muw = Aw.'*mu.';

% check that whitened sigma is the identity matrix
Sigmaw = Aw.'*Sigma*Aw; 
vertices = [0 1 1 0; 0 0 1 1];
transformed_vertices = Aw.'*vertices;

% find the new polygon's area (which is the new 1/pdf of the uniform dist)
A = polyarea(transformed_vertices(1,:),transformed_vertices(2,:));

% plot everything
x1 = -4:.1:4; x2 = -2.2:.1:6;
[X1, X2] = meshgrid(x1,x2);
Fw2 = mvnpdf([X1(:), X2(:)], muw.', Sigmaw);
Fw2 = reshape(Fw2, length(x2), length(x1));
figure()
contour(x1,x2,Fw2,[.0001 .001 1/A]);
line([transformed_vertices(1,:) transformed_vertices(1,1)],...
    [transformed_vertices(2,:), transformed_vertices(2,1)],...
    'linestyle','--','color','k');
text(muw(1), muw(2), '\omega_{2}'); text(muw(1) + 2, muw(2),'\omega_{1}');
text(2, 3.5, '\omega_{2}');
xlabel('x_{1}'); ylabel('x_{2}');
axis equal
hold off
figure()
surf(x1,x2,Fw2);
hold on
v = [transformed_vertices.', 1/A*ones(4,1)];
f = 1:4;
patch('Faces',f,'Vertices',v,'EdgeColor','b',...
    'FaceColor','none','linewidth',2);
hold off
xlabel('x_{1}'); ylabel('x_{2}'); zlabel('Probability Density');

%% Q3
%-------------------------------------------------------------------------%
%                      optimal desicion boundary
%-------------------------------------------------------------------------%
clear all %#ok
mu1 = zeros(2,1); mu2 = [0.75 0.5].';
SIGMA = [0.65 0.35;
    0.35 0.65];
R1 = mvnrnd(mu1, SIGMA, 1000);
R2 = mvnrnd(mu2, SIGMA, 1000);

x = linspace(-2,2,100);
y = linspace(-2,2,100);
[X, Y] = meshgrid(x,y);
Z = [];
for ii = 1:100
    z1 = mvnpdf([x(ii)*ones(100,1) Y(:,ii)],mu1.',SIGMA);
    z2 = mvnpdf([x(ii)*ones(100,1) Y(:,ii)],mu2.',SIGMA);
    Z = [Z (z1-z2)]; %#ok
end
figure()
contour(X,Y,Z)
hold on

% Find optimal Decision Boundary (using the formulas derived in the lecture
% notes)
w = SIGMA^-1*(mu1-mu2);
x0 = 0.5*(mu1+mu2);
n = cross([0 0 1].',[w;0]);
n = n(1:2)/norm(n);
p1 = x0-n*3; p2 = x0+n*3;
plot([p1(1),p2(1)],[p1(2),p2(2)],'--m','Linewidth',2)
xlim([-2 2]);ylim([-2 2]);
axis equal

%-------------------------------------------------------------------------%
%               optimal desicion boundary based on estimation
%-------------------------------------------------------------------------%

% Distributions estimation
muhat1 = 1/length(R1)*sum(R1);
muhat2 = 1/length(R2)*sum(R2);
murep1 = repmat(muhat1,1000,1);
murep2 = repmat(muhat2,1000,1);

% Assuming R1 and R2 have the same length, covariance matrix estimation
Sigmahat1 = (R1-murep1).'*(R1-murep1)/999;
Sigmahat2 = (R2-murep2).'*(R2-murep2)/999;

% Optimal Decision Based on Estimates
Sigmahat = (Sigmahat1+Sigmahat2)/2;
% sigmahat1 ~= sigmahat2 but close, therefore take the mean (we are asked
% to assume that they are equal)
w = Sigmahat^-1*(muhat1.'-muhat2.');
x0 = 0.5*(muhat1.'+muhat2.');
n = cross([0 0 1].',[w;0]);
n = n(1:2)/norm(n);
p1 = x0-n*3; p2 = x0+n*3;
plot([p1(1),p2(1)],[p1(2),p2(2)],'--r','Linewidth',2)
xlim([-2 2]);ylim([-2 2]);
legend('Probability Difference','Optimal Decision Boundary',...
    'Estimated Decision Boundary','Location','Southwest')
hold off
figure()
plot(R1(:,1),R1(:,2),'r^',R2(:,1),R2(:,2),'bv');
hold on

%-------------------------------------------------------------------------%
%               classification using statistical approach                  
%-------------------------------------------------------------------------%
[X, Y] = meshgrid(linspace(-4,4,100), linspace(-4,4,100));
X = X(:); Y = Y(:);
group = [ones(length(R1),1); 2*ones(length(R2),1)];
[C, err, ~, logp, coeff] = classify([X, Y], [R1; R2], group, 'linear'); 
K = coeff(1,2).const;
L = coeff(1,2).linear;
% plot stuff
gscatter(X,Y,C,'rb','.',3,'off');
f = @(x,y) K + [x y]*L;
h2 = ezplot(f,[-4 4 -4 4]);
set(h2,'Color','m','LineWidth',2);
axis([-4 4 -4 4]);
title('');
xlabel('x_{1}'); ylabel('x_{2}');
hold off

%-------------------------------------------------------------------------%
%                       classification using NN
%-------------------------------------------------------------------------%

P = [R1;R2].';
T = [ones(1000,1); zeros(1000,1)].';
net = newp([0 1; -2 2],1,'hardlim','learnwh');
P_Err = [];
lr = linspace(0.1,1e-5,50);
for ii = 1:50
    net.inputWeights{1}.learnParam.lr = lr(ii);
    net.Biases{1}.learnParam.lr = lr(ii);
    [net,Y,E] = adapt(net,P,T);
    P_Err = [P_Err, E*E.']; %#ok
end
figure();
plot(R1(:,1), R1(:,2), 'r^', R2(:,1), R2(:,2), 'bv');
legendstr = {'C1', 'C2'};
legend(legendstr);
plotpc(net.IW{1},net.b{1});
xlabel('x_1')
ylabel('x_2')
figure()
plot(P_Err,'--o')
ylabel('Sum Square Error')
xlabel('Epoch')

%-------------------------------------------------------------------------%
%                               verification
%-------------------------------------------------------------------------%

R1 = mvnrnd(mu1, SIGMA, 2000);
R2 = mvnrnd(mu2, SIGMA, 2000);
T = [ones(2000,1); -ones(2000,1)];
% because the covariance matrices are equal, the boundary is linear;
% therefore, in order to classify, we can use sign(v(x)) where v(x) =
% dot(w,x) + bias (just like a regular perceptron)

w_perceptron      = net.IW{1};          % Q3.d
b_perceptron      = net.b{1};
w_classifyFunc    = L.';                % Q3.c
b_classifyFunc    = K;
w_classifyOptimal = (w./norm(w)).';     % Q3.b
b_classifyOptimal = norm(x0);

y_perceptron      = sign([R1; R2]*w_perceptron.' + repmat(b_perceptron,4000,1));
y_classifyFunc    = sign([R1; R2]*w_classifyFunc.' + repmat(b_classifyFunc,4000,1));
y_classifyOptimal = sign([R1; R2]*w_classifyOptimal.' + repmat(b_classifyOptimal,4000,1));

sse_perceptron      = (T - y_perceptron).'*(T - y_perceptron);
sse_classifyFunc    = (T - y_classifyFunc).'*(T - y_classifyFunc);
sse_classifyOptimal = (T - y_classifyOptimal).'*(T - y_classifyOptimal);

disp(['SSE of Peceptron is ' num2str(sse_perceptron)]);
disp(['SSE of Classify function is ' num2str(sse_classifyFunc)]);
disp(['SSE of Optimal classification formula is ' num2str(sse_classifyOptimal)]);

%% Q4 Bayesian Classification On Half-Rings
[A, B] = RandInRing(10,6,1,1000);
%-------------------------------------------------------------------------%
%               classification using statistical approach                  
%-------------------------------------------------------------------------%
[X, Y] = meshgrid(linspace(-20,30,1000), linspace(-20,20,1000));
X = X(:); Y = Y(:);
group = [ones(length(A),1); 2*ones(length(B),1)];
[~, ~, ~, ~, coeff] = classify([X, Y], [A; B], group, 'linear'); 
w0_B = coeff(1,2).const;
w_B = coeff(1,2).linear;

%-------------------------------------------------------------------------%
%                       classification using NN
%-------------------------------------------------------------------------%
X = [ones(1000,1),A;
    ones(1000,1),B].';
T = [ones(1000,1);-ones(1000,1)].';
w_n = 1.5*rand(3,1)-ones(3,1);
Jp = Inf;
e = zeros(2000,1).';
eta = 1e-3;
while Jp(end) > 0
    ni = w_n.'*X;
    u = sign(ni);
    e = T-u;
    dw = eta*e*X.';
    Jp = [Jp, -0.5*ni*e.']; %#ok
    w_n = w_n + dw.';
end

%-------------------------------------------------------------------------%
%                       Testing the Classifiers
%-------------------------------------------------------------------------%
[A, B] = RandInRing(10,6,1,1000);
X = [A;B].';
T = [ones(1000,1);-ones(1000,1)].';
u_B = sign(w_B.'*X + w0_B);
E_B = (T-u_B)*(T-u_B).';
u_n = sign(w_n(2:3).'*X + w_n(1));
E_n = (T-u_n)*(T-u_n).';

% Plot Stuff
figure()
e_B = T-u_B;
C1good = find(e_B(1:1000) == 0);
plot(A(C1good,1),A(C1good,2),'bo')
hold on
FAv = find(e_B(1:1000) == 2);
plot(A(FAv,1),A(FAv,2),'ms')
C2good = find(e_B(1001:2000) == 0);
plot(B(C2good,1),B(C2good,2),'rx')
MDv = find(e_B(1001:2000) == -2);
plot(B(MDv,1),B(MDv,2),'g*')
plotpc(w_B.',w0_B(1))
axis([-20 30 -20 20])
axis equal
grid on
xlabel('x_1')
ylabel('x_2')
legend('C_1','False Alarm','C_2','Misdetection','Decision Boundary')
hold off