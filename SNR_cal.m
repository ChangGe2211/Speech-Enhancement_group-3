function snr=SNR_cal(clean,processed)

len=size(processed,1);   
clean=clean(1:len);   
Ps=sum(sum((clean-mean(mean(clean))).^2));
Pn=sum(sum((clean-processed).^2));            
snr=10*log10(Ps/Pn);



