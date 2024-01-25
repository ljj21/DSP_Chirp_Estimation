# DSP大作业

## 利用STFT验证信号参数与要求相同

MATLAB代码如下：

```matlab
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
```

使用了`spectrogram`函数进行了STFT，结果如下：

<img src="C:\Users\liuji\Desktop\DSP大作业\DSP大作业\1.png" alt="1" style="zoom:50%;" />

可以发现，在0-5$\mu s$ ，信号的频率从15MHz增加到25MHz，在5-10$\mu s$，信号的频率从25MHz下降到15MHz。由于采样率只有50MHz，故最多只能分辨25MHz的频率，对于高于25MHz的频率，将会发生混叠。25-35MHz的频率会被混叠到-25~-15MHz，也即25~15MHz，与图像相符。

## 分析滑动窗长度对时频域分辨率的影响

下面是仿真代码：

```matlab
load('wave_data.mat');
chirp_signal2 = chirp(t, 20e6, T, 32e6, "linear", 90);
observed_signal = chirp_signal + chirp_signal2 + wave_data;
figure;
window_size = 128; % 窗口大小
overlap = round(0.95 * window_size); % 重叠大小
nfft = window_size; % FFT点数
spectrogram(observed_signal, window_size, overlap, nfft, fs, 'yaxis');
title("Chirp Signal (w=" + num2str(window_size) + ")");
```

仿真时，调整`window_size`的大小，可以生成不同时频分辨率的信号：

<img src="C:\Users\liuji\Desktop\DSP大作业\DSP大作业\2_32.png" alt="2_32" style="zoom: 33%;" /><img src="C:\Users\liuji\Desktop\DSP大作业\DSP大作业\2_48.png" alt="2_48" style="zoom:33%;" />

<img src="C:\Users\liuji\Desktop\DSP大作业\DSP大作业\2_64.png" alt="2_64" style="zoom:33%;" /><img src="C:\Users\liuji\Desktop\DSP大作业\DSP大作业\2_72.png" alt="2_72" style="zoom:33%;" />

<img src="C:\Users\liuji\Desktop\DSP大作业\DSP大作业\2_96.png" alt="2_96" style="zoom:33%;" /><img src="C:\Users\liuji\Desktop\DSP大作业\DSP大作业\2_128.png" alt="2_128" style="zoom:33%;" />

较短的滑动窗长度意味着更高的时域分辨率，可以更好地分辨信号中快速变化的瞬时特性。然而，较短的窗口长度将导致频谱的主瓣展宽，频域分辨率下降。

对于`wave_data`，它的频率大概是**线性函数加上一个正弦函数**，如果窗长过大，快变的正弦信号就会“混叠”，难以分辨；而如果窗长过短，过低的频率分辨率会使两个chirp信号难以分辨。

综上，窗长为48时较为合适。

## 带噪单chirp信号估计

**参考了TIM上的"Polynomial Chirplet Transform With Application to Instantaneous Frequency Estimation"论文**，主要是想使用chirplet transform（CT）来分析chirp信号

对信号$s(t)$的CT的定义如下：
$$
CT_s(t_0,\omega,\alpha;\sigma)=\int_{-\infty}^{+\infty}\overline{z}(t)w_\sigma(t-t_0)\exp(-j\omega t)dt
$$
其中，
$$
z(t) = s(t)+j\mathcal{H}\{s(t)\}\\
\overline{z}(t)=z(t)\Phi_\alpha^R(t)\Phi_\alpha^M(t,t_0)\\
\Phi_\alpha^R(t) =\exp(-j\alpha t^2/2)\\
\Phi_\alpha^M(t,t_0)=\exp(j\alpha t_0 t)\\
w_\sigma(t) =\dfrac{1}{\sqrt{2\pi}\sigma}\exp\left(-\dfrac{t^2}{2\sigma^2}\right)
$$
$z(t)$是$s(t)$的解析信号，$w_\sigma(t)$是高斯窗函数，$CT$其实是对$\overline{z}(t)w_\sigma(t-t_0)$的傅里叶变换，只需要考察$\overline{z}(t)w_\sigma(t-t_0)$的性质，就能推测出CT的性质。

+ 在 time–frequency distributions (TFD）上，$\Phi_\alpha^R(t)$ 是一个旋转（**R**otate）操作，如果原信号的TFD是直线的，那么乘上$\Phi_\alpha^R(t)$后，直线的斜率会减小$\alpha$。
+ $\Phi_\alpha^M(t,t_0)$ 是一个频移操作，它将 $\omega$ 处的频率分量平移到$\omega+\alpha t_0$频率处。
+ 注意窗函数

下面的图（**摘自原论文**）较为形象地展示了chirplet变换对TFD的操作，假设原来的chirp信号的频率为$\omega_0+\lambda_0t$ (图中的青色线)

![Polynomial_Chirplet_Transform_With_Application_to_Instantaneous_Frequency_Estimation](C:\Users\liuji\Desktop\Polynomial_Chirplet_Transform_With_Application_to_Instantaneous_Frequency_Estimation.png)

CT先将其旋转，变成绿线，注意到直线的斜率变成了$\lambda_0-\alpha$；再进行平移操作，变成红线。在带宽为$\sigma$的高斯窗内，青色线的带宽为$\lambda_0 \sigma+1/\sigma$；绿色线和红色线的带宽为$\sigma|\lambda_0-\alpha|+1/\sigma$。当$\alpha =\lambda_0$ 时，带宽会取得最小值，$1/\sigma$。这意味着CT能使得TFD的**能量最为集中**。所以，$|CT_s(t_0,\omega,\alpha;\sigma)|$在$(\omega,\alpha)=(\omega_0,\lambda_0)$ 处取得全局最大值。这是我们估计的理论基础。



