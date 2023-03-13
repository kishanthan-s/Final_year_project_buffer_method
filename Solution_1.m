%% Generating a common waveform for NCO operartion

% Convolution with the ps filter and plotting the waveform 
fs = 7.68;
ps_out = conv(tmwaveform_5MHz,ps5_7_68);
figure 
plot(linspace(-fs/2,fs/2,length(tmwaveform_5MHz)), 20*log10(abs(fftshift(fft(tmwaveform_5MHz)))))
hold on
plot(linspace(-fs/2,fs/2,length(ps_out)), 20*log10(abs(fftshift(fft(ps_out)))))

% Interpolation
fs = 491.52;
tmw_interp491 = interp(tmwaveform_5MHz, 64);
figure 
plot(linspace(-fs/2,fs/2,length(tmw_interp491)), 20*log10(abs(fftshift(fft(tmw_interp491)))))

% waveBeforeNCO = conv(tmw_interp491,ps5_7_68);
% figure 
% plot(linspace(-fs/2,fs/2,length(waveBeforeNCO)), 20*log10(abs(fftshift(fft(waveBeforeNCO)))))

%% NCO operation for four bandwith parts with 5MHz bandwidth each
% center_frequency = 20; %MHz
% for ii = 1:length(ps_out)
%    NCO_out(ii) = tmw_interp491(ii)*exp(1i*2*pi*center_frequency/fs*ii);
% end

% Center frequencies for the bandwidth parts
center_freq1 = 20 ;%MHz
center_freq2 = 40 ;%MHz
center_freq3 = 60 ;%MHz
center_freq4 = 80 ;%MHz

% NCO operation is done seperately for each bandwidth part.
% for the ease of coding one for loop is used
for ii = 1:length(ps_out)
   NCO_out_1stBWP(ii) = tmw_interp491(ii)*exp(1i*2*pi*center_freq1/fs*ii);
   NCO_out_2ndBWP(ii) = tmw_interp491(ii)*exp(1i*2*pi*center_freq2/fs*ii);
   NCO_out_3rdBWP(ii) = tmw_interp491(ii)*exp(1i*2*pi*center_freq3/fs*ii);
   NCO_out_4thBWP(ii) = tmw_interp491(ii)*exp(1i*2*pi*center_freq4/fs*ii);
end

% Adding four BWPs to one spectrum
NCO_out = NCO_out_1stBWP + NCO_out_2ndBWP + NCO_out_3rdBWP + NCO_out_4thBWP;
figure
plot(linspace(-fs/2,fs/2,length(NCO_out)),20*log10(abs(fftshift(fft(NCO_out)))))

% Plotting the bandwith parts seperately
hold on;
plot(linspace(-fs/2,fs/2,length(NCO_out_1stBWP)),20*log10(abs(fftshift(fft(NCO_out_1stBWP)))))
hold on;
plot(linspace(-fs/2,fs/2,length(NCO_out_2ndBWP)),20*log10(abs(fftshift(fft(NCO_out_2ndBWP)))))
hold on;
plot(linspace(-fs/2,fs/2,length(NCO_out_3rdBWP)),20*log10(abs(fftshift(fft(NCO_out_3rdBWP)))))
hold on;
plot(linspace(-fs/2,fs/2,length(NCO_out_4thBWP)),20*log10(abs(fftshift(fft(NCO_out_4thBWP)))))


%% uplink operation
data = NCO_out;  % input data
temp_data = data; 

figure
plot(linspace(-fs/2,fs/2,length(temp_data)),20*log10(abs(fftshift(fft(temp_data)))))

%% NCO operation
% NCO Operation for 1st BWP
for ii = 1:length(ps_out)
    NCO_up_1stBWP(ii) =  temp_data(ii) * exp(-1i*2*pi*center_freq1/fs *ii);
end
figure
plot(linspace(-fs/2,fs/2,length(NCO_up_1stBWP)),20*log10(abs(fftshift(fft(NCO_up_1stBWP)))))
ps_out_uplink_BW1 = conv(ps5_7_68_uplink, NCO_up_1stBWP);
plot(linspace(-fs/2,fs/2,length(ps_out_uplink_BW1)), 20*log10(abs(fftshift(fft(ps_out_uplink_BW1)))))
% 
% NCO Operation for 2st BWP
for ii = 1:length(ps_out)
    NCO_up_2stBWP(ii) = temp_data(ii) * exp(-1i*2*pi*center_freq2/fs *ii);
end
figure
plot(linspace(-fs/2,fs/2,length(NCO_up_2stBWP)),20*log10(abs(fftshift(fft(NCO_up_2stBWP)))))
ps_out_uplink_BW2 = conv(ps5_7_68_uplink, NCO_up_2stBWP);
plot(linspace(-fs/2,fs/2,length(ps_out_uplink_BW2)), 20*log10(abs(fftshift(fft(ps_out_uplink_BW2)))))


% NCO Operation for 3st BWP
for ii = 1:length(ps_out)
    NCO_up_3stBWP(ii) = temp_data(ii) * exp(-1i*2*pi*center_freq3/fs *ii);
end
figure
plot(linspace(-fs/2,fs/2,length(NCO_up_3stBWP)),20*log10(abs(fftshift(fft(NCO_up_3stBWP)))))
ps_out_uplink_BW3 = conv(ps5_7_68_uplink, NCO_up_3stBWP);
plot(linspace(-fs/2,fs/2,length(ps_out_uplink_BW3)), 20*log10(abs(fftshift(fft(ps_out_uplink_BW3)))))


% NCO Operation for 4st BWP
for ii = 1:length(ps_out)
    NCO_up_4stBWP(ii) = temp_data(ii) * exp(-1i*2*pi*center_freq4/fs *ii);
end
figure
plot(linspace(-fs/2,fs/2,length(NCO_up_4stBWP)),20*log10(abs(fftshift(fft(NCO_up_4stBWP)))))
ps_out_uplink_BW4 = conv(ps5_7_68_uplink, NCO_up_4stBWP);
plot(linspace(-fs/2,fs/2,length(ps_out_uplink_BW4)), 20*log10(abs(fftshift(fft(ps_out_uplink_BW4)))))


 
 fs = 7.68;
 DECI_waveform_BW1 = decimate(ps_out_uplink_BW1,64);
 figure
 plot(linspace(-fs/2,fs/2,length(DECI_waveform_BW1)),20*log10(abs(fftshift(fft( DECI_waveform_BW1)))))
 title('Retrieved Waveform')
 xlabel('Frequency(MHz)')
 ylabel('Amplitude(dB)')
 
 DECI_waveform_BW2 = decimate(ps_out_uplink_BW2,64);
 figure
 plot(linspace(-fs/2,fs/2,length(DECI_waveform_BW2)),20*log10(abs(fftshift(fft( DECI_waveform_BW2)))))

 DECI_waveform_BW3 = decimate(ps_out_uplink_BW1,64);
 figure
 plot(linspace(-fs/2,fs/2,length(DECI_waveform_BW3)),20*log10(abs(fftshift(fft( DECI_waveform_BW3)))))

 DECI_waveform_BW4 = decimate(ps_out_uplink_BW1,64);
 figure
 plot(linspace(-fs/2,fs/2,length(DECI_waveform_BW4)),20*log10(abs(fftshift(fft( DECI_waveform_BW4)))))

 %checking the waveform whether they are equal
 if isequal(tmwaveform_5MHz, DECI_waveform_BW1)
    disp('The two waveforms are equal');
else
    disp('The two waveforms are not equal');
end
 clear temp_data; 
 
 [r, lag] = xcorr(tmwaveform_5MHz, DECI_waveform_BW1);
%plot the correlation
plot(lag, r);
title('Correlation vs Lag')
xlabel('Lag');
ylabel('Correlation');


%% 
x_time = ifft(NCO_out);
t_time = linspace(0,1,length(x_time));

figure
plot(t_time,real(x_time)); % Use "real" to remove imaginary part
title('Time Domain Waveform');
xlabel('Time (s)');
ylabel('Amplitude');
%%
% Example signal with cyclic prefix
signal_with_cp = [1 2 3 4 5 6 7 8 9 10];
cyclic_prefix_length = 3;
signal_length = length(signal_with_cp);

% Remove cyclic prefix
signal_without_cp = signal_with_cp(cyclic_prefix_length+1:signal_length);
% Display results
disp('Original signal with cyclic prefix:');
disp(signal_with_cp);
disp('Signal without cyclic prefix:');
disp(signal_without_cp);


