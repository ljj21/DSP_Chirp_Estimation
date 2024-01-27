+ `Chirplet_Transform.m`为chirplet变换的代码（**直接采用现有实现**（https://ww2.mathworks.cn/matlabcentral/fileexchange/72303-chirplet-transform）。
+ `extract_chirp_component.m`主要涉及对缩放去噪后的时频图做Hough变换，以及根据结果拟合直线斜率，合并预测结果。
+ `main.m`是主程序，可以分段运行，得到实验报告中每一个部分的结果。
+ `multi_chirp.m`主要涉及合并Hough变换后的结果。
+ `single_chirp.m`是用于对单chirp信号估计的函数。
+ `truncated_signal.m`是将信号在时域上截断。