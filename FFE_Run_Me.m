% About: This script plots the FFE coefficients to equalize the channel
% Author: Tejasvi Anand
% Date: 11/19/2015

clear all
%close all

%ChanFileName = 'TEC_Whisper42p8in_Meg6_THRU_C8C9.s4p';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the channel impulse response
% Load the S parameters
%Backplane = sparameters(ChanFileName);
%Data = Backplane.Parameters;
%Freq = Backplane.Frequencies;


load RationalFunc.mat RationalFunc
%load RationalFunc_RG58_Channel1.mat RationalFunc

npoles = length(RationalFunc.A); %returns largest array dimension
fprintf('The derived rational function contains %d poles.\n', npoles);
%Resp = freqresp(RationalFunc,Freq);% frequency response
%% Plot the frequency response
%figure;
%semilogx(Freq, 10.*log10(abs(Resp)));
%grid on;

%% Create a pulse response
Drate = 20e9; % Data rate = 5Gb/s 
Period = 1/Drate;
Trisefall = 10e-12;    % Rise/fall time [ps]
Tsample = 1e-12;     % sampling time 1ps 
Tsim = 20e-9;         % simulation time
Nsample = round(1/Drate/Tsample);
Nrisefall = round(Trisefall/Tsample);
InputPulse_one = [ (0:1:Nrisefall - 1)/Nrisefall ...
                            ones(1, Nsample - 2*Nrisefall) ...
                            (Nrisefall - 1:-1:0)/Nrisefall ];          
InputSignal_one = [InputPulse_one zeros(1, round(Tsim/Tsample))];
[Tresp_one, time_one] = timeresp(RationalFunc, InputSignal_one, Tsample); %ploting affect of channel on the impulse response?%


% Number of points in time
Points = round(Tsim/Tsample);
time = 0:Tsample:Tsim+Period-Tsample;

set(0,'DefaultTextFontSize',30) % Sets font size

%figure;
%plot(time.*1e12, InputSignal_one, 'b', 'LineWidth',2)
%hold on;
%plot(time.*1e12, Tresp_one, 'r', 'LineWidth',2)


%% Generate the PRBS 15 data

% Initial string for PRBS 15 data
init = [1 0 0 1 0 0 1 0 0 0 0 0 0 0 0];
g = [15 14];
z = prbs(init,g);

%% Create a waveform with the PRBS data
data_out = zeros(1,length(z)*Nsample + length(Tresp_one));
data_out = data_out';


for i=1:length(z)

    if z(i)==1
        data_out(1+(i-1)*Nsample:length(Tresp_one)+(i-1)*Nsample) = ...
            data_out(1 + (i-1)*Nsample : length(Tresp_one)+(i-1)*Nsample) + ...
            Tresp_one;
    else
        data_out(1 + (i-1)*Nsample : length(Tresp_one)+(i-1)*Nsample) = ...
            data_out(1 + (i-1)*Nsample : length(Tresp_one)+(i-1)*Nsample) - ...
            Tresp_one;
    end
end

%figure;
%plot(data_out, 'b', 'LineWidth',2)

%% Clip the raw data to plot the eye diagram
data_out_clip = zeros(1,(length(z)-1)*Nsample);
% Calcluate the index of the maximum value of impulse response
% This helps to estimate as how much delayed is the pulse response

[maxVal maxIndex] = max(Tresp_one);
data_out_clip = data_out(maxIndex : maxIndex+(length(z)-1)*Nsample);

%figure;
%plot(data_out_clip, 'b', 'LineWidth',2)

%h = eyediagram(data_out_clip,Nsample) %
%% FFE coefficients are estimated here

% Number of FFE coefficients
FFE_coeff = 3;

% Number of coefficients for the channel pulse response
Chnl_coeff = 20;
CHNL = zeros(1,Chnl_coeff);

% Coffecients are pre, main, post1, post2
% Estimate channel coefficient
[maxVal maxIndex] = max(Tresp_one);
for i=1:Chnl_coeff
    CHNL(i) = Tresp_one(maxIndex - Nsample + Nsample*(i-1));
end
CHNL = CHNL';

% Create the convulution matrix
CHNL_conv_matrix = convmtx(CHNL,FFE_coeff);

% Estimate the output Y
Y = zeros(1,Chnl_coeff+FFE_coeff-1)';
Y(2) = 1;
%Y(3) = 1; %for 3 tap
FFE_num = CHNL_conv_matrix'*Y;
FFE_den = (CHNL_conv_matrix'*CHNL_conv_matrix);
FFE_den_inv = (FFE_den)^-1;
FFE = FFE_den_inv*FFE_num;
FFE_mag = sum(abs(FFE));
FFE = FFE./FFE_mag

%% N TAP FFE implemented here

% FFE coefficients
%FFE = [-0.01 0.7 -0.3];

Tresp_ffe = zeros(1,length(FFE)*Nsample + length(Tresp_one));
Tresp_ffe = Tresp_ffe';

for i=1:length(FFE)
    Tresp_ffe(1+(i-1)*Nsample:length(Tresp_one)+(i-1)*Nsample) = ...
        Tresp_ffe(1+(i-1)*Nsample:length(Tresp_one)+(i-1)*Nsample) + ...
        FFE(i)*Tresp_one;
end

figure;
plot(Tresp_one, 'b', 'LineWidth',2)
hold on;
plot(Tresp_ffe, 'r', 'LineWidth',2)

% Create data with new response

data_ffe = zeros(1,length(z)*Nsample + length(Tresp_ffe));
data_ffe = data_ffe';

for i=1:length(z)

    if z(i)==1
        data_ffe(1+(i-1)*Nsample:length(Tresp_ffe)+(i-1)*Nsample) = ...
            data_ffe(1 + (i-1)*Nsample : length(Tresp_ffe)+(i-1)*Nsample) + ...
            Tresp_ffe;
    else
        data_ffe(1 + (i-1)*Nsample : length(Tresp_ffe)+(i-1)*Nsample) = ...
            data_ffe(1 + (i-1)*Nsample : length(Tresp_ffe)+(i-1)*Nsample) - ...
            Tresp_ffe;
    end
end

%figure;
%plot(data_ffe, 'b', 'LineWidth',2)

%% Clip the raw data to plot the eye diagram
data_ffe_clip = zeros(1,(length(z)-1)*Nsample);
% Calcluate the index of the maximum value of impulse response
% This helps to estimate as how much delayed is the pulse response

[maxVal maxIndex] = max(Tresp_ffe);
data_ffe_clip = data_ffe(maxIndex : maxIndex+(length(z)-1)*Nsample);

figure;
%plot(data_ffe_clip, 'b', 'LineWidth',2)

%eyediagram(data_ffe_clip,Nsample)

%Eye diagram from the communication tool bar
%Create an eye diagram scope object

h = commscope.eyediagram('SamplingFrequency', 1/Tsample, ...
                         'SamplesPerSymbol',Nsample, ...
                        'SymbolsPerTrace', 2, ...
                         'MinimumAmplitude', -0.2, ...
                         'MaximumAmplitude', 0.2, ...
                         'AmplitudeResolution', 0.00100, ...
                         'MeasurementDelay', 0, ...
                         'PlotType', '2D Color', ...
                         'PlotTimeOffset', 0, ...
                         'PlotPDFRange', [0 1], ...
                         'ColorScale', 'linear', ...
                         'RefreshPlot', 'on');

%Update the eye diagram
update(h, data_ffe_clip);
grid on;
%cmap = jet(64);
%cmap(1,:) = [0,0,0];
analyze(h);
h.Measurements



