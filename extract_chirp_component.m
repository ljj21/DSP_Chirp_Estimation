function [start_time, init_freq, duration, bandwidth, alpha] = extract_chirp_component(image, fs, fLevel, peak_num)
    % 读取二值化的边缘图像
    % bwImage = imread('image.png');
    % bwImage = imbinarize(image); % 可选步骤，将图像二值化
    bwImage = image;
    
    % 执行霍夫变换
    [H, theta, rho] = hough(bwImage);

    % 在霍夫空间中找到峰值
    peaks = houghpeaks(H, peak_num); % 选择peak_num个峰值点
    
    % 检测直线
    lines = houghlines(bwImage, theta, rho, peaks);
    
    detected_line = zeros(length(lines), 4);
    detected_line_num = 0;
    
    for index = 1:length(lines)
        is_detected = 1;
        point1 = lines(index).point1;
        point2 = lines(index).point2;
        for index2 = 1:detected_line_num
            point3 = detected_line(index2, 1:2);
            point4 = detected_line(index2, 3:4);
            if (sum(abs(point1 - point3)) < 20) && (sum(abs(point2 - point4)) < 20)
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
    time_per_x = 1 / fs;
    freq_per_y = fs / fLevel; 
    % point1 = lines.point1;
    % point2 = lines.point2;
    start_time = detected_line(:, 1) * time_per_x;
    init_freq = zeros(detected_line_num, 1);
    % duration = zeros(detected_line_num, 1);
    % bandwidth = zeros(detected_line_num, 1);
    alpha = zeros(detected_line_num, 1);
    for k1 = 1:detected_line_num
        left = detected_line(k1, 1);
        right = detected_line(k1, 3);
        x = left:right;
        y = zeros(1, length(x));
        for k2 = detected_line(k1, 1):detected_line(k1, 3)
            range = round((k2-left)/(right-left)*(detected_line(k1, 4)-detected_line(k1, 2))+detected_line(k1, 2));
            top = min(range+5, fLevel/2);
            low = max(1, range-5);
            signal_amp = image(low:top, k2);
            [~, max_index] = max(signal_amp);
            y(k2-left+1) = max_index + low - 1;
        end
        init_freq(k1) = y(1) * freq_per_y;
        p = polyfit(x, y, 1);
        alpha(k1) = p(1) * freq_per_y / time_per_x;

    end
    % init_freq = detected_line(:, 2) * freq_per_y;
    duration = (detected_line(:, 3) - detected_line(:, 1)) * time_per_x;
    bandwidth = abs(alpha .* duration);
    if (peak_num == 1)
        start_time = start_time(1);
        init_freq = init_freq(1);
        duration = duration(1);
        bandwidth = bandwidth(1);
        alpha = alpha(1);
    end
    % bandwidth = (detected_line(:, 4) - detected_line(:, 2)) * freq_per_y;
    % alpha = bandwidth ./ duration;
    % bandwidth = abs(bandwidth);
    % 绘制检测到的直线
    % figure;
    % imshow(bwImage);
    % hold on;
    % for k = 1:length(lines)
    %     xy = [lines(k).point1; lines(k).point2];
    %     plot(xy(:, 1), xy(:, 2), 'LineWidth', 2, 'Color', 'green');
    % end
    % figure;
    % imshow(bwImage);
    % hold on;
    % for k = 1:detected_line_num
    %     xy = [detected_line(k, 1:2); detected_line(k, 3:4)];
    %     plot(xy(:, 1), xy(:, 2), 'LineWidth', 2, 'Color', 'green');
    % end
    % title('Detected Lines');