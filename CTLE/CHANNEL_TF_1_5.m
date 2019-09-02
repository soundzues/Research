%Description: Fixed T/F of channel
%Problem: fix gain
%Fixed

clear all
close all

ChanFileName = 'TEC_Whisper42p8in_Meg6_THRU_C8C9.s4p';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the channel impulse response
% Load the S parameters
Backplane = sparameters(ChanFileName);
Data = Backplane.Parameters;
Freq = Backplane.Frequencies;
Z0 = Backplane.Impedance;


load RationalFunc_32.mat RationalFunc

npoles = length(RationalFunc.A);

%Get numerator and denominator for tf
A = RationalFunc.A;
C = RationalFunc.C;
den = cell(size(A));
num = cell(size(A));
k = 1;                          % Index of poles and residues
n = 0;                          % Index of numerators and denominators
while k <= npoles
    if isreal(A(k))             % Real poles
        n = n + 1;
        num{n} = C(k);
        den{n} = [1, -A(k)];
        k = k + 1;
    else                        % Complex poles
        n = n + 1;
        real_a = real(A(k));
        imag_a = imag(A(k));
        real_c = real(C(k));
        imag_c = imag(C(k));
        num{n} = [2*real_c, -2*(real_a*real_c+imag_a*imag_c)];
        den{n} = [1, -2*real_a, real_a^2+imag_a^2];
        k = k + 2;
    end
end
den = den(1:n);
num = num(1:n);

%trial -> chanel tf = 1/poles Q) Should I keep orignal num or num = 1?
%num = 1;

sys_channel = tf(num, den) %create transfer function

k = 1;
ntf = length(sys_channel);
H = sys_channel(1);

while k < ntf
    k = k + 1
    H = H*sys_channel(k)
    %H.num
end

%H.num
%h1 = sys_channel(1)

%h2 = sys_channel(2)

%h_tot = h1*h2


%{

%Get poles and zeros from RationalFunc
P = RationalFunc.A; %Pole
C = RationalFunc.C; %Zero

%}

%Create Zero-Pole-gain model from extracted RationalFunc

%sys_channel = zpk(1,P',1)



%Method 1: extracting Zeros Poles to create tf
%[z,p,k] = zpkdata(sys_channel)
%tf_channel = tf(1,p)

%Method 2: Converting zpk model to tf model
%[zero,pole] = zp2tf(0,P',1)
%tf_channel = tf(zero, Pole)

%Method 3: Model Reducer, tries to approixmate channel by rdeducing number
%of poles used
%modelReducer(sys_channel)

%resp = freqresp(sys_channel, Freq);

%figure(3)
%plot(Freq/1e9,20*log10(abs(resp)),'r');


%start at unity gain
%setting s to 0
gain = evalfr(H, 0);

%divide tf by gain 
H = H/gain;

%bode(h_tot)
%}


%Plot Poles and Zeros plot
figure(1)
pzmap(H)
grid on


%Plot bode plot
figure(2)
%Method 4: Use bode command
bode(H)
grid on

%Method 5: freqs
%freqs(0,P,Freq)

%Method 6: Semilog
%subplot(2,1,1);
%semilogx(Freq, 20.*log10(abs(sys_channel)));


%grid on
%}