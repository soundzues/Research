% Description: bode plots and use RationalFunc directly to create tf, start
% at unity gain with 100 poles
% Problem: Start at unity gain
% Fix (In Progress): Fix unity gain for CTLE and Channel
%                    Use ZPK to create CTLE TF

clear all
close all

%{
ChanFileName = 'TEC_Whisper42p8in_Meg6_THRU_C8C9.s4p';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the channel impulse response
% Load the S parameters
Backplane = sparameters(ChanFileName);
Data = Backplane.Parameters;
Freq = Backplane.Frequencies;
Z0 = Backplane.Impedance;

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
  'NPoles', [500 900]); %creates transfer function c/(s-a)
npoles = length(RationalFunc.A);
%}

load RationalFunc.mat RationalFunc


%Get numerator and denominator for tf
A = RationalFunc.A;
C = RationalFunc.C;

%trial -> chanel tf = 1/poles Q) Should I keep orignal num or num = 1?
%num = 1;

%create transfer function for channel
sys_channel = zpk(A', C', 1) 

%start at unity gain
%[Z,gain] = zero(sys_channel)

%setting s to 0
gain = evalfr(sys_channel, 0);

%divide tf by gain 
sys_channel = sys_channel/gain;

%Plot Poles and Zeros plot
figure(1)
pzmap(sys_channel)
grid on

%Plot bode plot
figure(2)
bode(sys_channel)
grid on


%descretprion: Place zero and pole atleast a decade apart for real values

%create transfer function for CTLE:
% H = (s+z)/((s+p1)(s+p2))
%put zero on location of poles 
zero_ctle = [1 5.732*10^7];

%chosing 2 poles at random locations > zero such that poles are at least 1
%decade apart from the zero
% p1 = 16*10^7
% p2 = 16.5*10^7

%pole_ctle = [1 1.45*10^8 5.25*10^15];

%pole_ctle = 1;

%pole_ctle = [1 250 15000]

%correct one:
%
pole_ctle = [1 3.25*10^8 2.64*10^16]


sys_ctle = tf(zero_ctle,pole_ctle)

%start at unity gain
%[Z,gain] = zero(sys_channel)

%setting s to 0
gain = evalfr(sys_ctle, 0);

%divide tf by gain 
sys_ctle = sys_ctle/gain;

%Plot Poles and Zeros plot
figure(3)
pzmap(sys_ctle)
grid on

%Plot Bode plot
figure(4)
%bode(sys_ctle)
margin(sys_ctle)
grid on


%CTLE and Channel poles and zeros:
figure(6)
pzmap(sys_channel, 'r', sys_ctle, 'g')
grid on

%create total response
sys_total = sys_channel*sys_ctle


%Plot Poles and Zeros plot
figure(7)
pzmap(sys_total)
grid on


%Bode plot
figure(8)
bode(sys_total)
grid on

%Plot frequency reponse
