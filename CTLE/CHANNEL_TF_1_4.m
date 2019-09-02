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

%{
% Convert to 2-port Diff. parameters
DiffData = s2sdd(Data);
DiffZ0 = 2*Z0;
% By default, |s2sdd| expects ports 1 & 3 to be inputs and ports 2 & 4 to
% be outputs. However if your data has ports 1 & 2 as inputs and ports 3 &
% 4 as outputs, then use 2 as the second input argument to |s2sdd| to
% specify this alternate port arrangement. For example,
% diffdata = s2sdd(data,2);

% Do the DC extrapolation if necessary
if(Freq(1) ~= 0)
  s21_dc = min(1, interp1(Freq, abs(squeeze(DiffData(2, 1, :))), 0, 'spline', 'extrap'));
  if(s21_dc == 1)
    s11_dc = 0;
  else
    s11_dc = max(0, min(sqrt(1 - (s21_dc.^2)), ...
      interp1(Freq, abs(squeeze(DiffData(1, 1, :))), 0, 'spline', 'extrap')));
  end
  
  DiffDataNew = cat(3, [s11_dc s21_dc; s21_dc s11_dc], DiffData);
  Freq = [0; Freq];
  % Compute the transfer function
  DiffTransFunc = s2tf(DiffDataNew,DiffZ0,DiffZ0,DiffZ0);
else
  % Compute the transfer function
  DiffTransFunc = s2tf(DiffData,DiffZ0,DiffZ0,DiffZ0);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Get the rational function fit and verify

RationalFunc = rationalfit(Freq,DiffTransFunc,-40, 'IterationLimit', 80, ...
  'NPoles', [0 100]); %creates transfer function c/(s-a)
npoles = length(RationalFunc.A);
%}

load RationalFunc_33.mat RationalFunc

%Get poles and zeros from RationalFunc
P = RationalFunc.A; %Pole
C = RationalFunc.C; %Zero


%Create Zero-Pole-gain model from extracted RationalFunc

sys_channel = zpk(1,P',1)

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
%gain = evalfr(sys_channel, 0)

%divide tf by gain 
%sys_channel = sys_channel/gain;



%Plot Poles and Zeros plot
figure(1)
pzmap(sys_channel)
grid on

%Plot bode plot
figure(2)

%Method 4: Use bode command
bode(sys_channel)

%Method 5: freqs
%freqs(0,P,Freq)

%Method 6: Semilog
%subplot(2,1,1);
%semilogx(Freq, 20.*log10(abs(sys_channel)));


grid on
