function [start_time, init_freq, duration, bandwidth, detected_alpha] = multi_chirp(observed_signal, fs, fLevel)
    % 对信号首先做一次chirplet变换，此时alpha = 0
    [Spec, ~, ~] = Chirplet_Transform(observed_signal, fLevel, 64, fs, 0);
    spec_abs = abs(Spec);
    max_spec_abs = max(max(spec_abs));
    % 去除噪声
    spec_abs(spec_abs < 0.5 * max_spec_abs) = 0;
    % 转化成图像，以便能利用图像处理中的技术
    scaled_data = uint8(255 * mat2gray(spec_abs));
    % 调用extract_chirp_component提取出不超过30条直线
    [~, ~, ~, ~, alpha2] = extract_chirp_component(scaled_data, fs, fLevel, 30);
    % 初始化估计参数的矩阵
    start_time = zeros(1, length(alpha2));
    init_freq = zeros(1, length(alpha2));
    duration = zeros(1, length(alpha2));
    bandwidth = zeros(1, length(alpha2));
    detected_alpha = zeros(1, length(alpha2));
    % 有效估计参数的个数
    valid_cnt = 0;
    % 遍历所有直线
    for k = 1:length(alpha2)
        isvalid = 1;
        % 调整alpha值为直线的斜率，使得对应信号的能量更集中
        alpha = alpha2(k);
        % 每个信号的迭代次数
        iter_cnt = 0;
        alpha_pre = alpha;
        % 当直线的斜率收敛时迭代停止
        while (iter_cnt == 0 || abs(alpha - alpha_pre) / abs(alpha_pre) > 0.02)
            alpha_pre = alpha;
            % 与函数开头的操作类似，只不过这里的alpha由0改成了直线的斜率
            [Spec, ~, ~] = Chirplet_Transform(observed_signal, fLevel, 256, fs, alpha);
            spec_abs = abs(Spec);
            max_spec_abs = max(max(spec_abs));
            spec_abs(spec_abs < 0.55 * max_spec_abs) = 0;
            scaled_data = uint8(255 * mat2gray(spec_abs));
            % 只提取一个直线
            [start_time_, init_freq_, duration_, bandwidth_, alpha] = extract_chirp_component(scaled_data, fs, fLevel, 1);
            iter_cnt = iter_cnt + 1;
            % 超过10次就不再迭代
            if(iter_cnt > 10)
                break;
            end
        end
        % 去除重复的迭代结果
        for index = 1:valid_cnt
            % 信号持续时间过短，或者估计的参数的误差都在10%以内就要删去
            if (duration_ < 2e-6 ...
                    || (abs(init_freq(index) - init_freq_) / abs(init_freq(index)) < 0.1 ...
                    && abs(duration(index) - duration_) / abs(duration(index)) < 0.1 ...
                    && abs(detected_alpha(index) - alpha) / abs(alpha) < 0.1))
                isvalid = 0;
                break;
            end

        end
        % 若新检测出来的参数和之前的没有重复，则更新矩阵
        if (isvalid == 1)
            valid_cnt = valid_cnt + 1;
            start_time(valid_cnt) = start_time_;
            init_freq(valid_cnt) = init_freq_;
            duration(valid_cnt) = duration_;
            bandwidth(valid_cnt) = bandwidth_;
            detected_alpha(valid_cnt) = alpha;
        end

    end
    % 把有效的估计量提取出来
    start_time = start_time(1:valid_cnt);
    init_freq = init_freq(1:valid_cnt);
    duration = duration(1:valid_cnt);
    bandwidth = bandwidth(1:valid_cnt);
    detected_alpha = detected_alpha(1:valid_cnt);
