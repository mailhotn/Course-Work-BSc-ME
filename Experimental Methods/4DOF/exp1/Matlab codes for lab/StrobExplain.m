Fs=41;
t=0:.0001:1;dt=t(2)-t(1);

y=sin(40*2*pi*t);
plot(t,y,'.-',t(1:round(1/Fs/dt):end),y(1:round(1/Fs/dt):end),'r')