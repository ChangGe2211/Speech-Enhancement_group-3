[y,Fs] = audioread('C:\Users\hpj19\Music\clean speech.wav');
[y2,Fs2] = audioread('C:\Users\hpj19\Music\stationary speech-shaped noise.wav' );

window=boxcar(100); %矩形窗
window1=hamming(100); %海明窗
window2=blackman(100); %blackman窗
noverlap=20; %数据无重叠
range='onesided'; %频率间隔为[0 Fs/2]，只计算一半的频率
nfft=1024;
[Pxx,f]=pwelch(10000*y(66506:66605),window,noverlap,nfft,Fs,range);
[Pxx1,f]=pwelch(10000*y(66506:66605),window1,noverlap,nfft,Fs,range);
[Pxx2,f]=pwelch(10000*y(66506:66605),window2,noverlap,nfft,Fs,range);
 
plot_Pxx=10*log10(Pxx);
plot_Pxx1=10*log10(Pxx1);
plot_Pxx2=10*log10(Pxx2);
 
figure(1)
plot(f/(Fs/2)*pi,plot_Pxx);
xlabel ('w/rad')
ylabel ('dB')
 
figure(2)
plot(f/(Fs/2)*pi,plot_Pxx1);
xlabel ('w/rad')
ylabel ('dB')
 
figure(3)
plot(f/(Fs/2)*pi,plot_Pxx2);  
xlabel ('w/rad')
ylabel ('dB')