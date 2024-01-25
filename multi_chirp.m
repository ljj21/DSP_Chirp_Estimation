function [start_time, init_freq, duration, bandwidth, detected_alpha] = multi_chirp(observed_signal, signal_num, fs, fLevel)

    [Spec, ~, ~] = Chirplet_Transform(observed_signal, fLevel, 64, fs, 0);
    spec_abs = abs(Spec);
    max_spec_abs = max(max(spec_abs));
    spec_abs(spec_abs < 0.5 * max_spec_abs) = 0;
    scaled_data = uint8(255 * mat2gray(spec_abs));
    [~, ~, ~, ~, alpha2] = extract_chirp_component(scaled_data, fs, fLevel, 10 * signal_num);
    start_time = zeros(1, length(alpha2));
    init_freq = zeros(1, length(alpha2));
    duration = zeros(1, length(alpha2));
    bandwidth = zeros(1, length(alpha2));
    detected_alpha = zeros(1, length(alpha2));
    valid_cnt = 0;

    for k = 1:length(alpha2)
        isvalid = 1;
        alpha = alpha2(k);
        iter_cnt = 0;
        alpha_pre = alpha;
        while (iter_cnt == 0 || abs(alpha - alpha_pre) / abs(alpha_pre) > 0.02)
            alpha_pre = alpha;
            [Spec, ~, ~] = Chirplet_Transform(observed_signal, fLevel, 256, fs, alpha);
            spec_abs = abs(Spec);
            max_spec_abs = max(max(spec_abs));
            spec_abs(spec_abs < 0.55 * max_spec_abs) = 0;
            scaled_data = uint8(255 * mat2gray(spec_abs));
            [start_time_, init_freq_, duration_, bandwidth_, alpha] = extract_chirp_component(scaled_data, fs, fLevel, 1);
            iter_cnt = iter_cnt + 1;
            if(iter_cnt > 10)
                break;
            end
        end

        for index = 1:valid_cnt

            if (duration_ < 2e-6 ...
                    || (abs(init_freq(index) - init_freq_) / abs(init_freq(index)) < 0.05 ...
                    && abs(duration(index) - duration_) / abs(duration(index)) < 0.05 ...
                    && abs(detected_alpha(index) - alpha) / abs(alpha) < 0.05))
                isvalid = 0;
                break;
            end

        end

        if (isvalid == 1)
            valid_cnt = valid_cnt + 1;
            start_time(valid_cnt) = start_time_;
            init_freq(valid_cnt) = init_freq_;
            duration(valid_cnt) = duration_;
            bandwidth(valid_cnt) = bandwidth_;
            detected_alpha(valid_cnt) = alpha;
        end

    end

    start_time = start_time(1:valid_cnt);
    init_freq = init_freq(1:valid_cnt);
    duration = duration(1:valid_cnt);
    bandwidth = bandwidth(1:valid_cnt);
    detected_alpha = detected_alpha(1:valid_cnt);
