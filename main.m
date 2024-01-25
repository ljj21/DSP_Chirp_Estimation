clear;
clc;
close all;
fs = 50e6; % 采样率为50MHz
T = 10e-6; % 脉冲持续时间为10μs
t = 0: 1/fs: T - 1/fs; % 时间向量
f0 = 15e6; % 起始频率为15MHz
f1 = 35e6; % 结束频率为35MHz

% 使用chirp函数生成线性调频信号
chirp_signal = chirp(t, f0, T, f1, 'linear', 90);
figure;
% 使用spectrogram函数进行短时傅里叶变换
window_size = 64; % 窗口大小
overlap = 32; % 重叠大小
nfft = window_size; % FFT点数
spectrogram(chirp_signal, window_size, overlap, nfft, fs, 'yaxis');
title('Chirp Signal (Spectrogram)');
%%
load('wave_data.mat');
chirp_signal2 = chirp(t, 20e6, T, 32e6, "linear", 90);
observed_signal = chirp_signal + chirp_signal2 + wave_data;
figure;
window_size = 72; % 窗口大小
overlap = round(0.95 * window_size); % 重叠大小
nfft = window_size; % FFT点数
spectrogram(observed_signal, window_size, overlap, nfft, fs, 'yaxis');
title("Chirp Signal (w=" + num2str(window_size) + ")");

%%
start_time = 3.8e-6;
duration = 7.9e-6 - start_time;
zero_freq = 6e6;
end_freq = 23.8e6;
chirp_signal = chirp(t, zero_freq, T, end_freq, "linear", 90);
init_freq = start_time / 1e-5 * (end_freq - zero_freq) + zero_freq;
bandwidth = abs(duration / 1e-5 * (end_freq - zero_freq));
observed_signal = truncated_signal(chirp_signal, start_time, start_time + duration, fs);
SNR = 3:30;
simu_num = 10;
start_time_est = zeros(simu_num, length(SNR));
init_freq_est = zeros(simu_num, length(SNR));
duration_est = zeros(simu_num, length(SNR));
bandwidth_est = zeros(simu_num, length(SNR));

for k = 1:length(SNR)
    for j = 1:simu_num
        noisy_signal = awgn(observed_signal, SNR(k));
        [start_time_est(j, k), init_freq_est(j, k), duration_est(j, k), bandwidth_est(j, k)] = single_chirp(noisy_signal, t, fs);
    end
end
% 计算RMSE
start_time_RMSE = sqrt(mean((start_time_est - start_time).^2, 1));
init_freq_RMSE = sqrt(mean((init_freq_est - init_freq).^2, 1));
duration_RMSE = sqrt(mean((duration_est - duration).^2, 1));
bandwidth_RMSE = sqrt(mean((bandwidth_est - bandwidth).^2, 1));
% 在一个图窗中绘制四个图像
figure;
subplot(2, 2, 1);
semilogy(SNR, start_time_RMSE);
xlabel('SNR/dB');
ylabel('RMSE');
title('Start Time RMSE');
subplot(2, 2, 2);
semilogy(SNR, init_freq_RMSE);
xlabel('SNR/dB');
ylabel('RMSE');
title('Initial Frequency RMSE');
subplot(2, 2, 3);
semilogy(SNR, duration_RMSE);
xlabel('SNR/dB');
ylabel('RMSE');
title('Duration RMSE');
subplot(2, 2, 4);
semilogy(SNR, bandwidth_RMSE);
xlabel('SNR/dB');
ylabel('RMSE');
title('Bandwidth RMSE');



%%

chirp_signal = chirp(t, 6e6, T, 23.8e6, "linear", 90);
chirp_signal2 = chirp(t, 9.7e6, T, 6e6, "linear", 90);
chirp_signal3 = chirp(t, 2e6, T, 15e6, "linear", 90);

observed_signal1 = truncated_signal(chirp_signal, 1e-6, 4e-6, fs);
observed_signal2 = truncated_signal(chirp_signal2, 4e-6, 7e-6, fs);
observed_signal3 = truncated_signal(chirp_signal3, 6e-6, 9e-6, fs);

observed_signal = observed_signal1 + observed_signal3 + observed_signal2;
% observed_signal = chirp(t, 4e6, T, 4e6, "linear", 90);
% observed_signal = awgn(observed_signal, 40);
signal_num = 3;
fLevel = 512;
[start_time, init_freq, duration, bandwidth, detected_alpha] = multi_chirp(observed_signal, 4, fs, fLevel);

for k = 1:length(start_time)
    fprintf('Start time of signal %d is %fe-6s\n', k, 1e6 * start_time(k));
    fprintf('Initial frequency of signal %d is %fe6Hz\n', k, 1e-6 * init_freq(k));
    fprintf('Duration of signal %d is %fe-6s\n', k, 1e6 * duration(k));
    fprintf('Bandwidth of signal %d is %fe6Hz\n', k, 1e-6 * bandwidth(k));
    fprintf('Alpha of signal %d is %fe11Hz/s\n', k, 1e-11 * detected_alpha(k));
    fprintf('\n');
end
