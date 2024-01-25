function [start_time, init_freq, duration, bandwidth] = single_chirp(observed_signal, t, fs)
    fLevel = 64;
    WinLen = 256;
    alpha = 0;
    alpha_pre = 0;
    iter_cnt = 0;

    while (iter_cnt == 0 || abs(alpha - alpha_pre) / abs(alpha_pre) > 0.02)
        alpha_pre = alpha;
        [Spec, Freq, ~] = Chirplet_Transform(observed_signal, fLevel, WinLen, fs, alpha);
        spec_abs = abs(Spec);
        [max_if_per_inst, index] = max(spec_abs);
        max_if = max(max_if_per_inst);
        observed_win = max_if_per_inst > 0.5 * max_if;
        observed_time = t(observed_win);
        observed_freq = Freq(index(observed_win));
        coefficients = polyfit(observed_time, observed_freq, 1);
        alpha = coefficients(1);
        iter_cnt = iter_cnt + 1;
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
    % plot(t, max_if_per_inst);
    % figure;
    % plot(t, Freq(index));
