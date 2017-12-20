function Strob()
A=1;DC=0;F=42;N=350;dstr=1.03;

% function Dan_Strob(F,A,N,DC,dstr)
% Generate sin signal while the Stroboscope is turned on 
% [A] = V;          Amplitude of input excitation;
% [DC] = V;         Additional DC to the signal;
% [F] = Hz;         Sine frequency;
% [N] = #;          Number of cycles
% [dstr]= #;        Strobe frequency deviation ~1.03

% Define Parameters
Fs = 5e3;             %  Sample frequency  Hz
Astr = 5;             % Amplitude of strobe signal
A = abs(A);
%% Create an NI session object and add two analog input channels  and two analog output channels
s = daq.createSession('ni');

%% daq.getDevices
% DevID = daqhwinfo('nidaq','InstalledBoardIds');    %% ID Number of NI USB device
% Dev = char(DevID);
d = daq.getDevices;
Dev = d.ID;

s.addAnalogInputChannel(Dev, 0:4, 'Voltage');
for ii=1:4, 
    s.Channels(ii).InputType='SingleEnded';
    %% Min Range is -0.20 to +0.20 Volts
    %% Max Range is -10   to 10 Volts
    s.Channels(ii).Range = [-5 5];
end 
s.Channels(5).InputType='SingleEnded';
%% Max Rate is 50 kHz in the current NI configuration.
s.Rate = Fs;

%% Output signal
% np = round(Fs./fn);
s.addAnalogOutputChannel(Dev, [0,1], 'Voltage');

%% Add listener for continuous data acquiring
lh = addlistener(s,'DataAvailable', @ ContPlotData);    % Plot continuous data

tt = 0:1/Fs:N/F; 
x = DC + A*sin(2*pi*F*tt.');
As = Astr.*ones(size(x));               % Strobe signal amplitude
xr(:,1) = x;
xr(:,2) = As+As.*square(2*pi*dstr*F*tt.',5); % strobe Signal 

s.queueOutputData(xr);         
s.NotifyWhenDataAvailableExceeds = Fs;	 % update graph every sec
%% Generate signal - Blocks Matlab while doing so
[~, ~] = s.startForeground;
delete(lh)

end
