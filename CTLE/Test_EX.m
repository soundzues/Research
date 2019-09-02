close all;
clear;

filename = 'connectedprobes.s2p';
orig_data = read(rfdata.data, filename);
freq = orig_data.Freq;
data = orig_data.S_Parameters(2,1,:);

%% rationalfit takes long -> do once, save the result in .mat, and load it

RationalFunc = rationalfit(freq, data,'NPoles',[0 10]);
save('RationalFunc_15in_12.mat', 'RationalFunc');

%% Compare the rational function with the original s2p file.

[resp, freq] = freqresp(RationalFunc, freq);
figure;
title('Rational fitting of S21 magnitude')
plot(orig_data,'S21','dB')
hold on
plot(freq/1e9,20*log10(abs(resp)),'r');