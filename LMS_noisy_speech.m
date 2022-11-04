
clear 
close all
clc
%% Generating noisy speech audios
[clean_speech_ori , fs] = audioread('clean speech.wav');
[sta_noise_ori, fs_sn] = audioread('stationary speech-shaped noise.wav');
[bab_noise_ori, fs_bn] = audioread('babble noise.wav');
SNR_1 = 1;
SNR_3 = 3;
SNR_5 = 5;
SNR_10 = 10;
N = 10*fs; % choosing the first 10 seconds of audio
clean_speech = clean_speech_ori(1:N);
sta_noise = sta_noise_ori(1:N);
bab_noise  = bab_noise_ori(1:N);
% clean_speech_1 = clean_speech - mean(clean_speech);
% clean_speech=clean_speech_1/max(abs(clean_speech_1));
noisy_speech_sta_1db = noisy_speech_generation(sta_noise,clean_speech,SNR_1,'noisy_speech_sta_1db.wav',fs);
noisy_speech_sta_3db = noisy_speech_generation(sta_noise,clean_speech,SNR_3,'noisy_speech_sta_3db.wav',fs);
noisy_speech_sta_5db = noisy_speech_generation(sta_noise,clean_speech,SNR_5,'noisy_speech_sta_5db.wav',fs);
noisy_speech_sta_10db = noisy_speech_generation(sta_noise,clean_speech,SNR_10,'noisy_speech_sta_10db.wav',fs);

noisy_speech_bab_1db = noisy_speech_generation(bab_noise,clean_speech,SNR_1,'noisy_speech_bab_1db.wav',fs);
noisy_speech_bab_3db = noisy_speech_generation(bab_noise,clean_speech,SNR_3,'noisy_speech_bab_3db.wav',fs);
noisy_speech_bab_5db = noisy_speech_generation(bab_noise,clean_speech,SNR_5,'noisy_speech_bab_5db.wav',fs);
noisy_speech_bab_10db = noisy_speech_generation(bab_noise,clean_speech,SNR_10,'noisy_speech_bab_10db.wav',fs);

t_win = 0.020; % 20ms window size
L_win = t_win*fs; 
hamming_win = hamming(L_win); % using hamming window

noisy_speech = noisy_speech_bab_3db;
L_noise = 10*L_win; % using the first 10 frames that only contains noise to esitimate noise 
noise_only = noisy_speech(1:L_noise);
% clean_speech_1 = clean_speech - mean(clean_speech);
% 
% clean_speech=clean_speech_1/max(abs(clean_speech_1));   
M=64;                                       
mu=0.001;  
itr=length(noisy_speech);                  
[y,W,e]=LMS(noisy_speech,clean_speech,M,mu,itr);
% filtered_speech=e/max(abs(e));
filtered_speech=e;
filter_name = 'LMS';
t=(0:N-1)/fs;
plot_results(filter_name, t,fs,clean_speech, noisy_speech,filtered_speech)
MSE =  MSE_cal(filtered_speech,clean_speech,fs);
stoi = stoi(clean_speech, filtered_speech, fs);
SNR_out=SNR_cal(clean_speech,filtered_speech);