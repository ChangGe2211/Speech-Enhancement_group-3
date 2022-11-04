clear 
close all
clc

%% This code for speech enhancement using wiener filter in frequency domain
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
% n_bart = floor(L_noise/L_win)-1;
% n_welch = floor(L_noise/ (L_win/ 2))- 1;
%bart(signal, L_win, window)
%%  spectrum estimation: Bartlett's method, Welch's method
P_noise_bart = bart(noise_only, L_win, hamming_win);
P_noise_welch = welch(noise_only,L_win, hamming_win);
 
% P_noisy_speech_bart = bart(noisy_speech, L_win, hamming_win);
% P_noisy_speech_welch = welch(noisy_speech,L_win, hamming_win);
% P_clean_speech_bart = bart(clean_speech, L_win, hamming_win);
% P_clean_speech_welch = welch(clean_speech, L_win, hamming_win);
N_frame_speech =  floor( length(noisy_speech)/ (L_win/ 2))- 1;
% N_frame_speech = floor(length(noisy_speech)/L_win);
out_fft = zeros(L_win,N_frame_speech);
n1 = 1;
for i = 1:N_frame_speech
    noise_speech_1 = noisy_speech(n1:n1+L_win-1).*hamming_win;
    clean_speech_1 = clean_speech(n1:n1+L_win-1).*hamming_win;
    fft_noise_speech = fft(noise_speech_1,L_win);
    fft_clean_speech = fft(clean_speech_1,L_win);
    P_noise_speech = abs(fft_noise_speech).^2;
    P_clean_speech = abs(fft_clean_speech).^2;
    Hw = P_clean_speech./(P_clean_speech +P_noise_welch);
    out_fft(:,i)=Hw.*fft_noise_speech;
    n1 = n1 +L_win/2;
end
filtered_speech = recover_signal(out_fft, fs);
MSE =  MSE_cal(filtered_speech,clean_speech,fs);
stoi = stoi(clean_speech, filtered_speech, fs);
SNR_out=SNR_cal(clean_speech,filtered_speech);
% Hw = P_clean_speech_bart./(P_clean_speech_bart +P_noise_bart);
% Hw = P_clean_speech_bart./(P_clean_speech_bart +P_noise_bart);
% P_filtered(:,i)=Hw.*y_fft;

t_min = 1*fs;
t_max = 3*fs;
min=-0.8;
max = 0.8;
filter_name = 'frequency domain weiner filter';
t=(0:N-1)/fs;
plot_results(filter_name, t,fs,clean_speech, noisy_speech,filtered_speech)

function plot_results(filter_name,t, fs,clean_speech, noisy_speech,filtered_speech)
% t=(0:N-1)/fs;
t_min = 1*fs;
t_max = 3*fs;
min=-0.8;
max = 0.8;
figure,
subplot(311);
plot(t(t_min:t_max),clean_speech(t_min:t_max));ylim([-max,max]);title('Clean speech');xlabel('t/s');ylabel('magnitude');
subplot(312);
plot(t(t_min:t_max),noisy_speech(t_min:t_max));ylim([-max,max]);title('Noisy speech');xlabel('t/s');ylabel('magnitude');
subplot(313);
plot(t(t_min:t_max),real(filtered_speech(t_min:t_max)));ylim([-max,max]);title(['Filtered speech: ',filter_name]);xlabel('t/s');ylabel('magnitude');
saveas(gcf, [filter_name,'audio.png']); 

figure,
subplot(311);
spectrogram(clean_speech(t_min:t_max),256,128,256,16000,'yaxis');title('Clean speech');
subplot(312);
spectrogram(noisy_speech(t_min:t_max),256,128,256,16000,'yaxis');title('Noisy speech');
subplot(313);
spectrogram(filtered_speech(t_min:t_max),256,128,256,16000,'yaxis');title(['Filtered speech: ',filter_name]);
saveas(gcf, [filter_name,'spec.png']); 
end

function P_bart = bart(signal0, L_win, window) %spectrum estimation bartlett's method
N_frame = floor(length(signal0)/L_win);
P_bart = 0;
n1 = 1;

for i = 1:N_frame
    signal = signal0(n1:n1+L_win-1).*window;
    P_peri = abs(fft(signal,L_win)).^2/length(signal);
    P_bart = P_bart + P_peri;
    n1 = n1 +L_win;
end
end

function P_welch = welch(signal0, L_win, window) %spectrum estimation bartlett's method
N_frame =  floor( length(signal0)/ (L_win/ 2))- 1;
P_welch = 0;
n1 = 1;
for i = 1:N_frame
    signal = signal0(n1:n1+L_win-1).*window;
    P_peri = abs(fft(signal,L_win)).^2/length(signal);
    P_welch = P_welch + P_peri;
    n1 = n1 +L_win/2;
end
end

%  
function signal = recover_signal(y, fs )

[freqRes, N_frame] = size(y);
win_len = freqRes;
inc = win_len / 2 ;
signal = zeros((N_frame - 1) * inc + win_len, 1);

for i = 1 : N_frame
    start = (i - 1) * inc + 1;
    spec = y(:, i);
    signal(start : start + win_len - 1) = signal(start : start + win_len - 1) + real(ifft(spec, win_len));
end
end
