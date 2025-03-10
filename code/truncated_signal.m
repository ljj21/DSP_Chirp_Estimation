function observed_signal = truncated_signal(signal, t1, t2, fs)
    % t1 为信号的起始时间
    % t2 为信号的结束时间
    % fs 为采样率
    len = length(signal);
    t = ones(1, len);
    index1 = floor(t1*fs) + 1;
    index2 = floor(t2*fs) + 1;
    t(1:index1 - 1) = 0;
    t(index2+1:end) = 0;
    observed_signal = signal .* t;

end