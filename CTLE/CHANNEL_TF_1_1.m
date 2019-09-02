%Description: analyze and fix chanel transfer function
%Prob: TF not working 


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
  'NPoles', [0 100]); %creates transfer function c/(s-a)
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


%Extract poles and zeros from rational function and place them in cell
%array
den = den(1:n);
num = num(1:n);

%convert cell array to ordinary array
%breaks cell array's into 3 seprate coloumn vectors 
num  = cell2mat(num) ; 
den = cell2mat(den) ;

%makes a single column vector of complex pole
num = complex(num(:,1) , num(:,2)) ; 
den = complex(den(:,2) , den(:,3)) ;

%
%den = [1 ; den] ;
%trial -> chanel tf = 1/poles Q) Should I keep orignal num or num = 1?
%num = 1;

den_tp = den';
num_tp = num';

%fprintf("test");
%sys_channel = tf(num_tp, den_tp) %create transfer function

sys_channel = tf(0,den_tp)

%P = pole(sys_channel);

%debug poles and zeros
%[num, den] = tfdata(sys_channel);
%[z,p,k] = tf2zp(num', den');

%Plot Poles and Zeros plot
figure(1)
pzmap(sys_channel)
grid on

%Plot bode plot
figure(2)
bode(sys_channel)
grid on
