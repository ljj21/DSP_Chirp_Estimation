function [start_time, init_freq, duration, bandwidth, alpha] = extract_chirp_component(image, fs, fLevel, peak_num)
    bwImage = image;
    % 执行霍夫变换，霍夫变换是图像处理中寻找直线的一种经典算法
    [H, theta, rho] = hough(bwImage);
    % 在霍夫空间中找到峰值
    peaks = houghpeaks(H, peak_num); % 选择peak_num个峰值点
    % 检测直线
    lines = houghlines(bwImage, theta, rho, peaks);
    % 初始化估计参数矩阵，每一行的格式为(x1,y1,x2,y2)
    detected_line = zeros(length(lines), 4);
    % 找到的直线很可能有重复，需要去重
    % 没有重复的直线数目
    detected_line_num = 0;
    % 遍历每一条Hough变换给出的直线
    for index = 1:length(lines)
        is_detected = 1;
        point1 = lines(index).point1;
        point2 = lines(index).point2;
        % 遍历之前找到的没有重复的直线
        for index2 = 1:detected_line_num
            point3 = detected_line(index2, 1:2);
            point4 = detected_line(index2, 3:4);
            % 如果两条直线的端点距离比较近，说明重复
            if (sum(abs(point1 - point3)) < 40) && (sum(abs(point2 - point4)) < 40)
                % 如果新直线的跨度比旧直线大，将旧直线更新为新直线
                if(point1(1) < point3(1) && point2(1) > point4(1))
                    detected_line(index2, :) = [point1, point2];
                end
                is_detected = 0;
                break;
            end
        end
        if(is_detected == 0)
            continue;
        end
        detected_line_num = detected_line_num + 1;
        detected_line(detected_line_num, :) = [point1, point2];
    end
    detected_line = detected_line(1:detected_line_num, :);
    % 横坐标一格代表多长时间
    time_per_x = 1 / fs;
    % 纵坐标一格代表多少频率
    freq_per_y = fs / fLevel; 
    start_time = detected_line(:, 1) * time_per_x;
    init_freq = zeros(detected_line_num, 1);
    alpha = zeros(detected_line_num, 1);
    % 接下来是拟合直线的斜率
    % 没有采用Hough变换给出的斜率，而是采用了与单chirp信号类似的最小二乘估计的方法
    for k1 = 1:detected_line_num
        % 把直线的左右端点找到
        left = detected_line(k1, 1);
        right = detected_line(k1, 3);
        x = left:right;
        y = zeros(1, length(x));
        % 扫描直线的每一个时刻
        for k2 = left:right
            % 对直线进行线性插值，确定纵坐标
            range = round((k2-left)/(right-left)*(detected_line(k1, 4)-detected_line(k1, 2))+detected_line(k1, 2));
            % 取一个带状区域[low:top]，防止取到其他信号的幅频响应
            % 预留了5的裕度
            top = min(range+5, fLevel/2);
            low = max(1, range-5);
            signal_amp = image(low:top, k2);
            [~, max_index] = max(signal_amp);
            y(k2-left+1) = max_index + low - 1;
        end
        init_freq(k1) = y(1) * freq_per_y;
        % 拟合斜率
        p = polyfit(x, y, 1);
        alpha(k1) = p(1) * freq_per_y / time_per_x;

    end
    duration = (detected_line(:, 3) - detected_line(:, 1)) * time_per_x;
    bandwidth = abs(alpha .* duration);
    % 只取一个峰
    if (peak_num == 1)
        start_time = start_time(1);
        init_freq = init_freq(1);
        duration = duration(1);
        bandwidth = bandwidth(1);
        alpha = alpha(1);
    end
    % 绘制检测到的直线
    % figure;
    % imshow(bwImage);
    % hold on;
    % for k = 1:detected_line_num
    %     xy = [detected_line(k, 1:2); detected_line(k, 3:4)];
    %     plot(xy(:, 1), xy(:, 2), 'LineWidth', 2, 'Color', 'green');
    % end
    % title('Detected Lines');