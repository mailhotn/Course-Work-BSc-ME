function poly3_fit(x,y)

% solution method: LS - using dtrend
yd = y;% detrend(y);
A = [0*x+1 x x.^2 x.^3];
q = A\yd; 

yfit = A*q; 
figure
h = plot(x,[yd yfit yd-yfit]);
set(h(1),'Marker','o','MarkerSize',6,'LineStyle','none')
set(h(2),'LineWidth',2)
set(h(3),'Marker','o','MarkerSize',6,'LineStyle','none','Color','r')
xlabel('$x$','interpreter','latex','FontSize',18)
ylabel('$y$','interpreter','latex','FontSize',18)
title(sprintf('$C_0=%.3f$ $C_1=%.3f$ $C_2=%.3f$ $C_3=%.3f$',q(1),q(2),q(3),q(4)),'interpreter','latex','FontSize',18)
hl = legend('$measured$','$fitted$','Error');
set(hl,'interpreter','latex','FontSize',18)
set(gcf,'Color','w')
set(gca,'FontSize',14)


end