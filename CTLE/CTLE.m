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
plot(data_out, 'b', 'LineWidth',2)

%{
%% Clip the raw data to plot the eye diagram
data_out_clip = zeros(1,(length(z)-1)*Nsample);
% Calcluate the index of the maximum value of impulse response
% This helps to estimate as how much delayed is the pulse response

[maxVal maxIndex] = max(Tresp_one);
data_out_clip = data_out(maxIndex : maxIndex+(length(z)-1)*Nsample);
%}

%figure;
%plot(data_out_clip, 'b', 'LineWidth',2)

%h = eyediagram(data_out_clip,Nsample) %


%%CTLE stuff 


%{
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
%}


