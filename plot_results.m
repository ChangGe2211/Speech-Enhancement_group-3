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
saveas(gcf, [filter_name,'_audio.png']); 

figure,
subplot(311);
spectrogram(clean_speech(t_min:t_max),256,128,256,16000,'yaxis');title('Clean speech');
subplot(312);
spectrogram(noisy_speech(t_min:t_max),256,128,256,16000,'yaxis');title('Noisy speech');
subplot(313);
spectrogram(filtered_speech(t_min:t_max),256,128,256,16000,'yaxis');title(['Filtered speech: ',filter_name]);
saveas(gcf, [filter_name,'_spec.png']); 
end

