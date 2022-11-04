function MSE = MSE_cal(filtered_speech0,clean_speech0,fs)
a = 0;

t_min = 1*fs;
t_max = 5*fs-1;

filtered_speech = filtered_speech0(t_min:t_max);
clean_speech = clean_speech0(t_min:t_max);
for it = 1:length(filtered_speech)
    Er  = filtered_speech (it) - clean_speech (it); % Error at point it
    a = a + Er^2;
    
end
MSE = sqrt(a/length(filtered_speech)); %Calculate minimum square error
end
% tt = 1:1: length(filtered);
% plot (Er.^2); 
% title ('Error');