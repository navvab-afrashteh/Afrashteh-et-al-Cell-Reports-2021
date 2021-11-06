function H = LPF_Gauss_25Hz(Fs)
if ~exist('Fs','var')
    Fs = 150;
end
if Fs == 150
    H = gaussmf(-4:4,[1.1 0]) ;
elseif Fs == 220
    H = gaussmf(-6:6,[1.6 0]) ;
end
H = H/sum(H);
