% Linearity script implements the linarity test as planned by the students
% in the preperatory work.
%
% The test tests the following linearity definition:
%
%  L[ a1*u1 + a2*u2] = a1*L[u1]+a2*L[u2]    (eq I) 
%
% we denote  L[ a1*u1 + a2*u2] as y_comb, L[u1] as y1
% and L[u2] as y2. So if eq I holds then the following 
% holds:
%
%      y_comb-a1*y1-a2*y2 = 0               (eq II)
%
% The script tests if eq II holds for the output of mass #1.

%% ~ insert here the parameters as desingned at home:
a1 = 5;
a2 = 4;
f1 = 5;
f2 = 3;

%% ~ define here the basic test parameters:
Fs = 5e3;            % Sampling frequency [Hz] (max. is 50e3 Hz)
dt = 1/Fs;           % Sampling time interval defined by the frequency
t0 = 0:dt:10;        % t0 is a time vector for the fisrt 10 seconds to start
                     % with a 10 seconds zero command 
t1 = (10+dt):dt:20;  % t1 is a time vector for the next 10 seconds for 
                     % the actual test.

t = [t0 t1];         % t is the total time vector for the test     

u0 = zeros(1,length(t0))';% u0 is the 10 seconds initial zero command
u1 = sin(2*pi*f1*t1)';    % u1 is a 10 seconds sin signal with frequency f1 
u2 = sin(2*pi*f2*t1)';    % u2 is a 10 seconds sin signal with frequency f2  

u1 = [u0 ;u1];        % add the zero command to first 10 seconds of u1
u2 = [u0 ;u2];        % add the zero command to first 10 seconds of u2
%% ~ Run the linearity test functions:
%  LinearityTest(Fs,user_defind_input,'File name')

LinearityTest(Fs,a1*u1+a2*u2,'CombOut'); % performs y_comb = L[ a1*u1 + a2*u2]
LinearityTest(Fs,u1,'y1Out');            % performs y1 = L[u1]
LinearityTest(Fs,u2,'y2Out');            % performs y2 = L[u2]

%% ~ Load the data from the tests and check eq II for mass #1:
load('CombOut');y_comb = data(:,1)-mean(data(:,1)); % loads data of mass #1 from y_comb
load('y1Out');y1 = data(:,1)-mean(data(:,1));       % loads data of mass #1 from y1
load('y2Out');y2 = data(:,1)-mean(data(:,1));       % loads data of mass #1 from y2
     
err = y_comb-(a1*y1+a2*y2);                         % if eq II holds then err~0


%% ~ Plot the results:

figure(1);

% plot the actual outputs and observe how they match
subplot 211
plot(time,y_comb,time,a1*y1+a2*y2)
xlabel 'Time [sec]'
ylabel 'Output of mass #1'
legend('y_{comb}','a1*y1+a2*y2')

% plot the error (eq II)
subplot 212
plot(time,err,'r')
xlabel 'Time [sec]'
ylabel 'Error'

