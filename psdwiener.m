[y2,Fs2] = audioread('stationary speech-shaped noise.wav');
[y,Fs] = audioread('clean speech.wav');
Length=length(y);
Length2=length(y2);
Start=66506;
SNR=5;
y2=y2/norm(y2,2).*10^(-SNR/20)*norm(y);
Voice=y(Start:Start+99)+y2(Start:Start+99);
y=y(Start:Start+99);
N=100;
Rxx=xcorr(Voice,N-1,'biased'); 
  for i=1:N
      for j=1:N
        mRxx(i,j)=Rxx(N-i+j); 
      end
  end
Rxd=xcorr(Voice,y,N-1,'biased'); 
for i=1:N
    mRxd(i)=Rxd(N-1+i); 
end
h = inv(mRxx)*mRxd'; % calculate wiener cofficients
a(1)=1;
for i=2:100
    a(i)=0;
end
b = h;
N = 2000;
x = randn(1,N);
y = filter(b,a,x,[],2);% true power spectrum
omega = 0:0.01:pi;
z = exp(1i*omega);
num = abs(polyval(b,z));
den = abs(polyval(a,z));
P_true = 1*num./den;
P_peri = abs(fft(y,1024)).^2/length(y); 
plot(linspace(0,pi,512),10*log10(P_peri(1:512))) %peridogram
hold on
plot(omega,10*log10(P_true),'--')% true power spectrum
xlabel ('w/rad')
ylabel ('dB')