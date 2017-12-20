function [ZETA] = damping_coeff(time,data,ff) 
% damping_coeff(t,X,NF)
%
% calculate the damping coefficient for a decaying sine.
% time - vector of smpled time (N) 
% data - matrix with 4 colums and N rows 
% ff - natural frequency [HZ]
%
% It is importent to choose good initial guess

plot(time,data(:,1:4))
title('choose 2 points on exp. decay part')
shg
[xx,yy]=ginput(2);

xx=sort(xx);
i1 = find(time>=xx(1) & time<=xx(2));

plot(time(i1),data(i1,1:4))
shg
t=time(i1);
figure;
for q1=1:4
        
    y = detrend(data(:,q1));
    y3 = y(i1);
%     plot(time,y,time(i1),y(i1))
    % fit decaying sine
%     drawnow

    ft = fittype('a+b*exp(c*x)*cos(d*x+e)');
    %coeffnames(ft)
    
    InitalGuess = [0,1/2,-1e-3*ff*2*pi,ff*2*pi,0.3];
    [crv,gof,output] = fit(t,y3,ft,'StartPoint',InitalGuess);

    yf = feval(crv,t);
    figure()
    plot(t,y3,t,yf ,t,output.residuals)
    legend('measured','fitted','error')

    title(sprintf('f_d=%.1f Hz, f_n=%.1f Hz, zeta=%.4f',crv.d/2/pi, sqrt(crv.c^2+crv.d^2)/2/pi,-crv.c/sqrt(crv.d^2+crv.c^2)))
    Fd(q1) = crv.d/2/pi;
    FN(q1) = sqrt(crv.c^2+crv.d^2)/2/pi;
    ZETA(q1) = -crv.c/sqrt(crv.d^2+crv.c^2);
end
disp('The damped N.F[Hz] are as follows:')
fprintf('%.3f\n',Fd.')
disp('The N.F[Hz] are as follows:')
fprintf('%.3f\n',FN.')
disp('The damping coefficients are as follows:')
fprintf('%.5f\n',ZETA.')
end