clear;
clc;
close all;
% load('wave_data.mat')
% 生成线性调频信号
fs = 50e6; % 采样率为50MHz
T = 10e-6; % 脉冲持续时间为10μs
t = 0:1 / fs:T - 1 / fs; % 时间向量

f0 = 15e6; % 起始频率为15MHz
f1 = 35e6; % 结束频率为35MHz

% 使用chirp函数生成线性调频信号
% chirp_signal = chirp(t, f0, T, f1, 'linear', 90); % 'linear'表示线性变化，90表示相位角度
% chirp_signal2 = chirp(t, 20e6, T, 32e6, "linear", 90);
% noisy_signal = awgn(chirp_signal2, 40);
% observed_signal = truncated_signal(chirp_signal2, noisy_signal,2e-6, 6e-6, fs);
% observed_signal = chirp_signal;
% 绘制线性调频信号的时域图
% figure;
% plot(t, real(chirp_signal));
% xlabel('Time (s)');
% ylabel('Amplitude');
% title('Chirp Signal (Time Domain)');

% 使用spectrogram函数进行短时傅里叶变换
% in part 1, window size is set to 64
window_size = 64; % 窗口大小
overlap = 63; % 重叠大小
nfft = window_size; % FFT点数

% figure;
% spectrogram(observed_signal, window_size, overlap, nfft, fs, 'yaxis');

% title('Chirp Signal (Spectrogram)');
%%

chirp_signal = chirp(t, 5e6, T, 25e6, "linear", 90);
chirp_signal2 = chirp(t, 4e6, T, 9.5e6, "linear", 90);
chirp_signal3 = chirp(t, 25e6, T, 15e6, "linear", 90);

observed_signal1 = truncated_signal(chirp_signal, 5e-6, 9e-6, fs);
observed_signal2 = truncated_signal(chirp_signal2, 3.65e-6, 6.3e-6, fs);
observed_signal3 = truncated_signal(chirp_signal3, 4.2e-6, 8.5e-6, fs);

observed_signal = observed_signal1 + observed_signal3 + observed_signal2;
% observed_signal = chirp(t, 4e6, T, 4e6, "linear", 90);
% observed_signal = awgn(observed_signal, 40);
signal_num = 3;
fLevel = 512;
[start_time, init_freq, duration, bandwidth, detected_alpha] = multi_chirp(observed_signal, 3, fs, fLevel);
