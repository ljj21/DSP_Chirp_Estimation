function [start_time, init_freq, duration, bandwidth] = single_chirp(observed_signal, t, fs)
    % 频率轴下标的个数
    fLevel = 256;
    % 高斯窗的窗长
    WinLen = 64;
    % 估计出来的频率变化率
    alpha = 0;
    alpha_pre = 0;
    % 迭代次数
    iter_cnt = 0;
    % 若前后两次估计出来的alpha相差太多，则继续迭代
    while (iter_cnt == 0 || abs(alpha - alpha_pre) / abs(alpha_pre) > 0.02)
        alpha_pre = alpha;
        % 进行chirplet变换
        % Spec为频率响应
        [Spec, Freq, ~] = Chirplet_Transform(observed_signal, fLevel, WinLen, fs, alpha);
        spec_abs = abs(Spec);
        % max_if_per_inst：每个时刻的幅频响应的最大值
        % index: 最大值对应的频率下标
        [max_if_per_inst, index] = max(spec_abs);
        % 整个时频图中的最大值
        max_if = max(max_if_per_inst);
        % 忽略最大幅频响应较低的时刻
        observed_win = max_if_per_inst > 0.5 * max_if;
        % 取出这些时刻，认为这就是chirp信号的持续时间
        observed_time = t(observed_win);
        % 取出这些时刻的幅频响应最大值对应的频率
        observed_freq = Freq(index(observed_win));
        % 使用线性最小二乘拟合
        coefficients = polyfit(observed_time, observed_freq, 1);
        % 斜率即为频率变化率
        alpha = coefficients(1);
        iter_cnt = iter_cnt + 1;
        % 至多迭代10次
        if(iter_cnt > 10)
            break;
        end
    end

    start_time = observed_time(1);
    init_freq = observed_freq(1);
    duration = observed_time(end) - observed_time(1);
    bandwidth = abs(observed_freq(end) - observed_freq(1));
    % fprintf('Start time: %d s\n', observed_time(1));
    % fprintf('Initial frequency: %d Hz\n', observed_freq(1));
    % fprintf('Duration: %d s\n', observed_time(end)-observed_time(1));
    % fprintf('Bandwidth: %d Hz\n', abs(observed_freq(end)-observed_freq(1)));
    

    % figure;
    % 绘制每个时刻的最大频响
    % plot(t, max_if_per_inst);
    % 绘制每个时刻最大频响对应的频率值
    % figure;
    % plot(t, Freq(index));
    % 绘制chirplet变换后的时频图
    % figure;
    % imagesc(t,Freq,abs(Spec));
    % JET = colormap(jet);
    % colormap (JET);   % Set the current colormap Style 
    % box on;
    % colorbar off;
