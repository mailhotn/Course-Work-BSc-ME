function damping_coeff(time,data,ff) 
% damping_coeff(t,X,NF)
%
% calculate the damping coefficient for a decaying sine.
% time - vector of smpled time (N) 
% data - matrix with 4 colums and N rows 
% ff - natural frequency [HZ]
%
% It is importent to choose good initial guess

figure(1)
plot(time,data(:,1:4))
title('choose 2 points on exp. decay part','interpreter','latex','FontSize',18)
set(gcf,'Color','w')
set(gca,'FontSize',14)
shg
[xx,~]=ginput(2);

xx = sort(xx);
i1 = find(time>=xx(1) & time<=xx(2));

% disply selected data
figure(2)
plot(time(i1),data(i1,1:4))
set(gcf,'Color','w')
set(gca,'FontSize',14)
shg
t = time(i1);

for q1=1:4,
        
    y = detrend(data(:,q1));        % remove DC and drift
    y3 = y(i1);
    
    figure(3)
    plot(time,y,time(i1),y(i1))
    set(gcf,'Color','w')
    set(gca,'FontSize',14)
    drawnow
    
    % fit decaying sine
    ft = fittype('a+b*exp(c*x)*cos(d*x+e)'); % Assumed model
    %coeffnames(ft)
    
    InitalGuess = [0,0.3,0,ff*2*pi,-pi];
    [crv,gof,output] = fit(t,y3,ft,'StartPoint',InitalGuess);

    yf = feval(crv,t);

    figure
    h = plot(t,[y3 yf output.residuals]);
    set(h,'LineWidth',2)
    set(h(1),'Marker','x')
    set(h(3),'Color','r')
    xlabel('time[sec]','interpreter','latex','FontSize',18)
    ylabel('$[V]$','interpreter','latex','FontSize',18)
    shg
    hl = legend('measured','fitted','error');
    set(hl,'interpreter','latex','FontSize',18)
    set(gcf,'Color','w')
    set(gca,'FontSize',14)
    title(sprintf('$f_d=%.3f$[Hz], $\\zeta=%.4f$',crv.d/2/pi,-crv.c/sqrt(crv.d^2+crv.c^2)),'interpreter','latex','FontSize',18)

    FN(q1) = crv.d/2/pi;
    ZETA(q1) =-crv.c/sqrt(crv.d^2+crv.c^2);
end
disp('The damped N.F[Hz] are as follows:')
fprintf('%.3f\n',FN)
disp('The damping coefficients are as follows:')
fprintf('%.5f\n',ZETA)
end