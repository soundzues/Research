%Description: bode plots and use RationalFunc directly to create tf (In progress)
% 2 complex poles system.

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
  'NPoles', [0 2]); %creates transfer function c/(s-a)
npoles = length(RationalFunc.A);

%2 extracted poles - conjugate 
% p1 = -5.732*10^7 + i*4.4675*10^9
% p2 = -5.732*10^7 - i*4.4675*10^9
% (s - p1)*(s - p2) = (s - (-5.732*10^7 + i*4.4675*10^9))*(s - (-5.732*10^7 - i*4.4675*10^9))
% s^2 + 1.1464*10^8 + 1.199618*10^19

%create transfer functions for channel

%zero = 1;

%poles = [1 1.1464*10^8 1.199618*10^19];

%poles = [RationalFunc.A]

%sys_channel = tf(zero, poles);

%sys_channel = tf(RationalFunc)


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


%{
  Channel Transfer function:
            1
  -------------------------
  s^2 + 1.146e08 s + 1.2e19
%}

%Plot Poles and Zeros plot
figure(1)
pzmap(sys_channel)
grid on

%Plot bode plot
figure(2)
bode(sys_channel)

%get freq response
resp = freqresp(sys_channel, Freq);

%create transfer function for CTLE:
% H = (s+z)/((s+p1)(s+p2))
%put zero on location of poles 
zero_ctle = [1 5.732*10^7];

%chosing 2 poles at random locations > zero
% p1 = 7.5*10^7
% p2 = 7.0*10^7
%(s+7.0*10^7)+(s+7.5*10^7) = (s^2 + 1.45*10^8s + 5.25*10^15)
pole_ctle = [1 1.45*10^8 5.25*10^15];

sys_ctle = tf(zero_ctle,pole_ctle);

%{
  CTLE Transfer function:
            1
  -------------------------
  s^2 + 1.146e08 s + 1.2e19
%}

%Plot Poles and Zeros plot
figure(4)
pzmap(sys_ctle)
grid on

%Plot Bode plot
figure(5)
bode(sys_ctle)

%CTLE and Channel poles and zeros:
figure(6)
pzmap(sys_channel, 'r', sys_ctle, 'g')
grid on

%create total response
sys_total = sys_channel*sys_ctle

%{
Equalized system Transfer function 
                 s + 5.732e07
  --------------------------------------------------------
  s^4 + 2.596e08 s^3 + 1.202e19 s^2 + 1.74e27 s + 6.298e34
%}


%Plot Poles and Zeros plot
figure(7)
pzmap(sys_total)
grid on


%Bode plot
figure(8)
bode(sys_total)


%Plot frequency reponse
