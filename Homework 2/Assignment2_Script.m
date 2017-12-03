%% Q2
mu = [0.5 0.5];
Sigma = [0.1 0.05;
    0.05 0.1];
x1 = -1:.03:2; x2 = x1;
[X1, X2] = meshgrid(x1,x2);
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
% find the new polygon's area 
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
clear all %#ok
mu1 = zeros(2,1); mu2 = [0.75 0.5].';
SIGMA = [0.65 0.35;
    0.35 0.65];
R1 = mvnrnd(mu1, SIGMA, 1000);
R2 = mvnrnd(mu2, SIGMA, 1000);

% Distributions estimation
muhat1 = 1/length(R1)*sum(R1);
muhat2 = 1/length(R2)*sum(R2);
Sigmahat1 = zeros(2);
Sigmahat2 = zeros(2);

% Assuming R1 and R2 have the same length, covariance matrix estimation
for ii = 1:length(R1)
    Sigmahat1 = Sigmahat1 + (R1(ii,:) - muhat1).'*(R1(ii,:) - muhat1);
    Sigmahat2 = Sigmahat2 + (R2(ii,:) - muhat2).'*(R2(ii,:) - muhat2);
end
Sigmahat1 = Sigmahat1/(length(R1)-1);
Sigmahat2 = Sigmahat2/(length(R2)-1);

plot(R1(:,1),R1(:,2),'r^',R2(:,1),R2(:,2),'bv');
hold on

%-------------------------------------------------------------------------%
%               classification using statistical approach
%-------------------------------------------------------------------------%
[X, Y] = meshgrid(linspace(-4,4,100), linspace(-4,4,100));
X = X(:); Y = Y(:);
group = [ones(length(R1),1); 2*ones(length(R2),1)];
[C, err, P, logp, coeff] = classify([X, Y], [R1; R2], group, 'linear'); 
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
[A, B] = RandInRing(10,6,1,1000);
P = [A; B].';
T = [ones(1000,1); zeros(1000,1)].';
figure()
plotpv(P,T)
net = perceptron;
linehandle = plotpc(net.IW{1},net.b{1});
E = 1;
perf = [];
ii = 1;
while (sse(E))
    ii = ii + 1;
    [net,Y,E,~,~,tr] = adapt(net,P,T);
    perf(1,ii) = tr.perf; %#ok
    linehandle = plotpc(net.IW{1},net.b{1},linehandle);
    drawnow;
end
figure();
plot(A(:,1), A(:,2), 'bo', B(:,1), B(:,2), 'rx');
legendstr = {'C1', 'C2'};
legend(legendstr);
plotpc(net.IW{1},net.b{1});
figure();
plot(1:ii, perf, '--o');
xlabel('Epoch'); ylabel('Performance (MSE)');
% Boundary validation
[vA, vB] = RandInRing(10,6,1,3000);
vP = [vA; vB].';
vT = [ones(3000,1); zeros(3000,1)].';
y = net(vP);
MD = 0; FA = 0;
MDind = []; FAind = [];
for ii = 1:length(y)
    switch y(ii) - vT(ii)
        case 1
            MD = MD + 1;
            MDind = [MDind; ii]; %#ok
        case -1
            FA = FA + 1;
            FAind = [FAind; ii]; %#ok
    end
end
figure()
plot(A(:,1), A(:,2), 'bo', B(:,1), B(:,2), 'rx');
hold on
legendstr = {'C1', 'C2'};
if ~isempty(FAind)
    plot(vP(1,FAind), vP(2,FAind),'ms', 'linewidth', 2);
    legendstr = [legendstr, 'False Alarm'];
end
if ~isempty(MDind)
    plot(vP(1,MDind), vP(2,MDind), 'gd', 'linewidth', 2)
    legendstr = [legendstr, 'Miss-detection'];
end
legend(legendstr);
plotpc(net.IW{1},net.b{1});
hold off