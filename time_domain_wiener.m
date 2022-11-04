[noi,Fs2] = audioread('C:\Users\hpj19\Music\babble noise.wav');
[y,Fs] = audioread('C:\Users\hpj19\Music\clean speech.wav');
clean = upsample(y,ceil(Fs2/Fs));
Length=length(clean);
Length2=length(noi);
Start=floor(rand*abs(Length2-Length));
Voice=clean;
SNR=1;       %Fixed signal noise ratio
noi=noi/norm(noi,2).*10^(-SNR/20)*norm(clean); %calculate noise with fixed SNR
Voice(1:Length)=clean+noi(Start+1:Start+Length); %Get the noisy signal
N=100; % Order of wiener filter

for k=1:5000
  Rxx=xcorr(Voice(N*(k-1)+1:N*k).*hamming(N),N-1,'biased');% choose window type
  for i=1:N
      for j=1:N
        mRxx(i,j)=Rxx(N-i+j); % N*N degree autocorrelation function
      end
  end
  Rxd=xcorr(Voice(N*(k-1)+1:N*k).*hamming(N),clean(N*(k-1)+1:N*k).*hamming(N),N-1,'biased'); 
  for i=1:N
    mRxd(i)=Rxd(N-1+i); % 1*N degree cross-correlation functions
  end
h = inv(mRxx)*mRxd'; % The optimal solution of the filter is obtained through wiener-Hopf formula, degree of h (N*1)
yy = filter(h,1,Voice(N*(k-1)+1:N*k)); 
filtered(N*(k-1)+1:N*k) = yy;  %Filter result
end
clean=clean(1:length(filtered));

plot_results('Time-domain wiener',(0:length(filtered)-1)/Fs, Fs,clean, Voice,filtered)
MSE_time_wiener = MSE_cal(filtered,clean,Fs);
SNR_out=SNR_cal(clean(1*Fs:3*Fs),transpose(filtered(1*Fs:3*Fs)));
STOI_out=stoi(clean(1*Fs:3*Fs),filtered(1*Fs:3*Fs),Fs);

% t=(0:length(Yy)-1)/Fs2;
% t1=(0:length(y3)-1)/Fs2;
% subplot(211)
% plot(t,Yy);
% xlabel ('t/s')
% ylabel ('filtered signal')
% subplot(212)
% plot(t1,y3);
% sound(Yy,Fs2)
figure(1)
N=1024;
P_peri = 1000*peri(clean(66506:66605),N);
plot(linspace(0,pi,N/2),10*log10(P_peri(1:N/2)))
xlabel ('w/rad')
ylabel ('dB')
figure(2)
N=1024;
P_peri = 1000*peri(Voice(66506:66605),N);
plot(linspace(0,pi,N/2),10*log10(P_peri(1:N/2)))
xlabel ('w/rad')
ylabel ('dB')
figure(3)
N=1024;
P_peri = 1000*peri(filtered(66506:66605),N);
plot(linspace(0,pi,N/2),10*log10(P_peri(1:N/2)))
xlabel ('w/rad')
ylabel ('dB')
% SP=1000*Yy;
% Nw = 32000; 							% Window function length
% window = blackman(32000);  			% Window function selection
% noverlap = 30000; 					% overlap length
% nfft = 2^nextpow2(length(window)); 	% DFT number
% fs = 16000; 
% spectrogram(SP(3*fs:30*fs), window, noverlap, nfft,fs, 'yaxis');
% subplot(212)
% plot(t,Voice)
