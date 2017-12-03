%% Q2
mu = [0.5,0.5].';
Sigma = [0.1 0.05;
         0.05 0.1];
PHI = 0.7071*[-1 1;
              1 1];
Lambda = 0.01*[5 0;
               0 15];
A = PHI*Lambda^(-0.5);

%% Q3
mu1 = [0 0].';
mu2 = [0.75 0.5].';
S = [0.65 0.35;
     0.35 0.65];
x = linspace(-2,2,100);
y = linspace(-2,2,100);
[X, Y] = meshgrid(x,y);
Z = [];
for ii = 1:100
    z1 = mvnpdf([x(ii)*ones(100,1) Y(:,ii)],mu1.',S);
    z2 = mvnpdf([x(ii)*ones(100,1) Y(:,ii)],mu2.',S);
    Z = [Z (z1-z2)];
end
contour(X,Y,Z)
hold on
% Find optimal Decision Boundary
w = S^-1*(mu1-mu2);
x0 = 0.5*(mu1+mu2);
n = cross([0 0 1].',[w;0]);
n = n(1:2)/norm(n);
p1 = x0-n*3; p2 = x0+n*3;
plot([p1(1),p2(1)],[p1(2),p2(2)],'--m','Linewidth',2)
xlim([-2 2]);ylim([-2 2]);
legend('Probability Difference','Optimal Decision Boundary','Location','Southwest')
axis equal
R1 = mvnrnd(mu1, S, 1000);
R2 = mvnrnd(mu2, S, 1000);
% Distributions estimation
muhat1 = 1/length(R1)*sum(R1);
muhat2 = 1/length(R2)*sum(R2);
murep1 = repmat(muhat1,1000,1);
murep2 = repmat(muhat2,1000,1);
Sigmahat1 = (R1-murep1).'*(R1-murep1)/999;
Sigmahat2 = (R2-murep2).'*(R2-murep2)/999;

Shat = (Sigmahat1+Sigmahat2)/2;
w = Shat^-1*(muhat1-muhat2);
x0 = 0.5*(muhat1+muhat2);
n = cross([0 0 1].',[w;0]);
n = n(1:2)/norm(n);
p1 = x0-n*3; p2 = x0+n*3;
plot([p1(1),p2(1)],[p1(2),p2(2)],'--r','Linewidth',2)
xlim([-2 2]);ylim([-2 2]);