close all;
clear;

%% Open channel s2p file, rationalfit it, or load it

%ChanFileName = 'Z:/Documents/MATLAB/Dicode_Math/Channel/15in_12.S2P';

%load('Z:/Documents/MATLAB/iPWM_Math/RationalFunc_15in_12.mat', 'RationalFunc');

ChanFileName = 'TEC_Whisper42p8in_Meg6_THRU_C8C9.s4p';
load RationalFunc_30.mat RationalFunc

orig_data = read(rfdata.data, ChanFileName);
freq = orig_data.Freq;
data = orig_data.S_Parameters(2,1,:);

%% rationalfit takes long -> do once, save the result in .mat, and load it

% RationalFunc = rationalfit(freq, data,'NPoles',[500 900]);
% save('RationalFunc_15in_12.mat', 'RationalFunc');

%% Compare the rational function with the original s2p file.

[resp, freq] = freqresp(RationalFunc, freq);
figure;
title('Rational fitting of S21 magnitude')
plot(orig_data,'S21','dB')
hold on
plot(freq/1e9,20*log10(abs(resp)),'r');