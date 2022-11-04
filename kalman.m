[clean,~]=audioread('C:\Users\hpj19\Music\clean speech.wav');      %read clean speech
[noise,fs]=audioread('C:\Users\hpj19\Music\babble noise.wav');       %read noise
N=5*fs;                                              
clean=clean(1:N); % only take first 5 seconds
noise=noise(1:N);
t=(0:N-1)/fs;
SNR=1;                 
noise=noise/norm(noise,2).*10^(-SNR/20)*norm(clean);
voice=clean+noise;          % add noise with fixed SNR
Time = (0:1/fs:(length(clean)-1)/fs)';      
Voice=voice(fs+1:2*fs,1);        % The first 1 second speech Voice used for variance estimate
len_win = 0.0015;       % length of window 1.5ms
shift_percent = 1;       % shift step
filter_order = 20;           % filter order
iter = 1;                      % iteration times

len_winframe = fix(len_win * fs);% adding window
window = ones(len_winframe,1);
HoppingSamples = fix(len_winframe.*shift_percent);

num_seg = fix(((length(voice)-len_winframe)/HoppingSamples) + 1);
Index = (repmat(1:len_winframe,num_seg,1) + repmat((0:(num_seg-1))'*HoppingSamples,1,len_winframe))';
y = voice(Index).*repmat(window,1,num_seg);

C = [zeros(1,filter_order-1),1];   % observation matrix
R = var(Voice);  
[filt_coeff, Q] = lpc(y, filter_order);              % LPC prediction getting coefficient of the filter 
P = R * eye(filter_order,filter_order);              % covariance matrix
filtered = zeros(1,length(voice));    % filtered signal
filtered(1:filter_order) = voice(1:filter_order,1)';   % initialize
x_old = voice(1:filter_order,1);

i = filter_order+1;
j = filter_order+1;
for k = 1:num_seg  
    begin = j;     
    OutputOld = x_old;    % Keep the first filter_order estimates for each iteration
    for l = 1:iter               % iteration times
        A = [zeros(filter_order-1,1) eye(filter_order-1); fliplr(-filt_coeff(k,2:end))];

        for l = i:len_winframe % kalman filter process
            x = A * x_old;
            P_old = (A * P * A') + (C' * Q(k) * C);
            K = (P_old * C')/((C * P_old * C') + R);
            x_old = x + (K * (y(l,k) - (C*x)));
            filtered(j-filter_order+1:j) = x_old';
            P = (eye(filter_order) - K * C) * P_old;
            j = j+1;
        end
        i = 1;
        if l < iter
            j = begin;
            x_old = OutputOld;
        end
        [filt_coeff(k,:), Q(k)] = lpc(filtered((k-1)*len_winframe+1:k*len_winframe),filter_order);    % lpc of updated signal
    end
end
plot_results('kalman',(0:length(filtered)-1)/fs, fs,clean, voice,filtered)
MSE_time_kalman = MSE_cal(filtered,clean,fs);
SNR_out=SNR_cal(clean(1*fs:3*fs),transpose(filtered(1*fs:3*fs)));
STOI_out=stoi(clean(1*fs:3*fs),filtered(1*fs:3*fs),fs);