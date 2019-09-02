%%% To generate and save the rational fit for s4p parameters
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
  'NPoles', [0 1]);
npoles = length(RationalFunc.A);
fprintf('The derived rational function contains %d poles.\n', npoles);
Resp = freqresp(RationalFunc,Freq);% frequency response

%% Part (b) Compare fitted result and S-parameter
figure;
subplot(2,1,1);
semilogx(Freq, 10.*log10(abs(DiffTransFunc)));
subplot(2,1,2);
semilogx(Freq, 10.*log10(abs(Resp)));
hold on; 
semilogx(Freq, angle(Resp)./pi.*180);
hold on;
semilogx(Freq, angle(DiffTransFunc)./pi.*180);

%% Part (c) Worst case eye diagram
Drate = 5e9;            % Data rate = 5Gb/s
Period = 1/Drate;
Trisefall = 20e-12;    % Rise/fall time [ps]
Tsample = 1e-12;     % sampling time 1ps 
Tsim = 15e-9;         % simulation time
Nsample = round(1/Drate/Tsample);
Nrisefall = round(Trisefall/Tsample);
InputPulse_one = [ (0:1:Nrisefall - 1)/Nrisefall ...
                            ones(1, Nsample - 2*Nrisefall) ...
                            (Nrisefall - 1:-1:0)/Nrisefall ];          
InputSignal_one = [InputPulse_one zeros(1, round(Tsim/Tsample))];
[Tresp_one, time_one] = timeresp(RationalFunc, InputSignal_one, Tsample);

% Number of points in time
Points = round(Tsim/Tsample);
time = 0:Tsample:Tsim+Period-Tsample;
%[ui, myeye] = worsteye(Tresp_one);
%negmyeye = -1.*myeye;


set(0,'DefaultTextFontSize',30) % Sets font size
figure;
plot(time.*1e12, InputSignal_one, 'b', 'LineWidth',2)
hold on;
plot(time.*1e12, Tresp_one, 'r', 'LineWidth',2)


%plot(ui, myeye);
%hold on;
%plot(ui, negmyeye);
%%

