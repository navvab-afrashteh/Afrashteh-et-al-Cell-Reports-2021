function Hd = BPF_LS(Fs,Fstop1,Fpass1,Fpass2,Fstop2,dur)
%BPF_LS Returns a discrete-time filter object.
% FIR least-squares Bandpass filter designed using the FIRLS function.

% All frequency values are in Hz.
if ~exist('Fs','var')
    Fs = 100;
elseif isempty(Fs)
    Fs = 100;
end
if ~exist('Fstop1','var')
    Fstop1 = 0.25;
elseif isempty(Fstop1)
    Fstop1 = 0.25;
end
if ~exist('Fpass1','var')
    Fpass1 = 0.5;
elseif isempty(Fpass1)
    Fpass1 = 0.5;
end
if ~exist('Fpass2','var')
    Fpass2 = 6;
elseif isempty(Fpass2)
    Fpass2 = 6;
end
if ~exist('Fstop2','var')
    Fstop2 = 6.5;
elseif isempty(Fstop2)
    Fstop2 = 6.5;
end
if ~exist('dur','var')
    dur = 5;
elseif isempty(dur)
    dur = 5;
end
N      = Fs*dur;   % Order
Wstop1 = 1;     % First Stopband Weight
Wpass  = 1;     % Passband Weight
Wstop2 = 1;     % Second Stopband Weight

% Calculate the coefficients using the FIRLS function.
b  = firls(N, [0 Fstop1 Fpass1 Fpass2 Fstop2 Fs/2]/(Fs/2), [0 0 1 1 0 ...
           0], [Wstop1 Wpass Wstop2]);
Hd = dfilt.dffir(b);

% [EOF]
