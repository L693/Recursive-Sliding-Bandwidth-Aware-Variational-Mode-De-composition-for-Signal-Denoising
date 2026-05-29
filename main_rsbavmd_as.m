
clear; clc; close all;


%% ====================== 1. 生成复合模拟信号 ======================

fs = 50;          % 采样频率 (Hz)
T = 60;            % 信号总时长 (s)
t = (0:1/fs:T-1/fs)'; % 时间向量（强制列向量）
N = length(t);     % 总采样点数

% 1.1 Damped vibration mode
s_d = 2* exp(-0.1 * t) .* cos(4 * pi * t + pi/4);

% 1.2 Intermodulation component
s_1 = 1.8 * cos(10 * pi * t) ;

s_2 = 0.8 * cos(20 * pi * t) ;

% % 1.3 Impulse noise
% s_i = zeros(size(t));
% t_impulse = 10; 
% impulse = 10 * exp(-50 * (t - t_impulse)) .* (t >= t_impulse);
% s_i = s_i + impulse;

% 合成最终信号
% x_clean = s_d + s_1  + s_2 + s_i; % 干净复合信号


x_clean = s_d + s_1  + s_2;

snr_input = 2;
x_noisy = awgn(x_clean, snr_input);

%% ====================== 2. 定义评价指标函数 ======================
% 信噪比计算
function snr_val = calculate_snr(clean, denoised)
    clean = clean(:);
    denoised = denoised(:);
    noise = clean - denoised;
    snr_val = 10*log10(mean(clean.^2)/mean(noise.^2));
end

% 均方误差计算
function mse_val = calculate_mse(clean, denoised)
    clean = clean(:);
    denoised = denoised(:);
    mse_val = mean((clean - denoised).^2);
end

% 均方根误差计算
function rmse_val = calculate_rmse(clean, denoised)
    clean = clean(:);
    denoised = denoised(:);
    rmse_val = sqrt(mean((clean - denoised).^2));
end

% 相关系数计算
function corr_val = calculate_corr(clean, denoised)
    clean = clean(:);
    denoised = denoised(:);
    corr_val = corrcoef(clean, denoised);
    corr_val = corr_val(1,2);
end

% 包络平滑度计算
function smooth_val = calculate_envelope_smoothness(signal)
    signal = signal(:);
    envelope = abs(hilbert(signal));
    smooth_val = -mean(abs(diff(envelope)));
end

% 排列熵计算
function pe = permutation_entropy(x, m, tau)
    x = x(:);
    n = length(x);
    
    permlist = perms(1:m);
    num_patterns = factorial(m);
    count = zeros(1, num_patterns);
    
    for i = 1:n-(m-1)*tau
        [~, idx] = sort(x(i:tau:i+(m-1)*tau));
        for j = 1:num_patterns
            if isequal(idx, permlist(j,:))
                count(j) = count(j) + 1;
                break;
            end
        end
    end
    
    count = count(count > 0);
    p = count / sum(count);
    pe = -sum(p .* log(p)) / log(num_patterns); % 归一化到[0,1]
end

% 多尺度排列熵计算
function mpe = multiscale_permutation_entropy(signal, max_scale, m, tau)
    signal = signal(:);
    mpe = zeros(1, max_scale);
    
    for scale = 1:max_scale
        % 粗粒化过程
        coarse_len = floor(length(signal)/scale);
        coarse_signal = zeros(coarse_len, 1);
        
        for i = 1:coarse_len
            coarse_signal(i) = mean(signal((i-1)*scale+1:i*scale));
        end
        
        % 计算该尺度下的排列熵
        mpe(scale) = permutation_entropy(coarse_signal, m, tau);
    end
end

%% ====================== 3. 实现所有对比算法 ======================

methods = {'Noisy', 'EMD', 'EEMD', 'Bandpass Filter', 'Grid-VMD', 'Wavelet', 'RSBAVMD'};
num_methods = length(methods);
results = zeros(num_methods, 6); % [SNR, MSE, RMSE, Correlation, Envelope Smoothness, Time(s)]
denoised_signals = cell(num_methods, 1);
iter_info = struct('vmd_iter', [], 'vmd_time', [], 'rsbavmd_iter', [], 'rsbavmd_time', []);

% 3.1 原始含噪信号
tic;
denoised_signals{1} = x_noisy;
time_noisy = toc;
results(1,:) = [calculate_snr(x_clean, x_noisy), calculate_mse(x_clean, x_noisy), ...
                calculate_rmse(x_clean, x_noisy), calculate_corr(x_clean, x_noisy), ...
                calculate_envelope_smoothness(x_noisy), time_noisy];

% 3.2 EMD算法
tic;
imf_emd = emd(x_noisy, 'Interpolation','pchip');
denoised_emd = sum(imf_emd(:,1:4), 2); % 选择前4个能量最大的IMF重构
time_emd = toc;
denoised_signals{2} = denoised_emd;
results(2,:) = [calculate_snr(x_clean, denoised_emd), calculate_mse(x_clean, denoised_emd), ...
                calculate_rmse(x_clean, denoised_emd), calculate_corr(x_clean, denoised_emd), ...
                calculate_envelope_smoothness(denoised_emd), time_emd];

% 3.3 EEMD算法（内嵌兼容版）
tic;
imf_eemd = eemd_compatible(x_noisy, 0.2, 100); % 噪声标准差0.2，集成次数100
denoised_eemd = sum(imf_eemd(:,1:4), 2);
time_eemd = toc;
denoised_signals{3} = denoised_eemd;
results(3,:) = [calculate_snr(x_clean, denoised_eemd), calculate_mse(x_clean, denoised_eemd), ...
                calculate_rmse(x_clean, denoised_eemd), calculate_corr(x_clean, denoised_eemd), ...
                calculate_envelope_smoothness(denoised_eemd), time_eemd];

% 3.4 理想带通滤波（已知有效频率范围0-20Hz）
tic;
[b,a] = butter(4, 20/(fs/2), 'low'); % 4阶巴特沃斯低通滤波器
denoised_bp = filtfilt(b,a,x_noisy); % 零相位滤波
time_bp = toc;
denoised_signals{4} = denoised_bp;
results(4,:) = [calculate_snr(x_clean, denoised_bp), calculate_mse(x_clean, denoised_bp), ...
                calculate_rmse(x_clean, denoised_bp), calculate_corr(x_clean, denoised_bp), ...
                calculate_envelope_smoothness(denoised_bp), time_bp];

% 3.5 网格优化VMD
K_range = 2:8;
alpha_range = 500:500:5000;
best_snr = -inf;
best_K = 4;
best_alpha = 2000;
for K = K_range
    for alpha = alpha_range
        [u,~,~,~] = VMD(x_noisy, alpha, 0, K, 0, 1, 1e-7);
        denoised = sum(u(1:K-1,:), 1)'; % 转置为列向量
        current_snr = calculate_snr(x_clean, denoised);
        if current_snr > best_snr
            best_snr = current_snr;
            best_K = K;
            best_alpha = alpha;
        end
    end
end
tic;
[u_vmd,~,~,vmd_iter] = VMD(x_noisy, best_alpha, 0, best_K, 0, 1, 1e-7);
denoised_vmd = sum(u_vmd(1:best_K-1,:), 1)';
time_vmd = toc;
denoised_signals{5} = denoised_vmd;
results(5,:) = [calculate_snr(x_clean, denoised_vmd), calculate_mse(x_clean, denoised_vmd), ...
                calculate_rmse(x_clean, denoised_vmd), calculate_corr(x_clean, denoised_vmd), ...
                calculate_envelope_smoothness(denoised_vmd), time_vmd];
iter_info.vmd_iter = vmd_iter;
iter_info.vmd_time = time_vmd;

% 3.6 小波去噪

tic;
[thr, sorh, keepapp] = ddencmp('den','wv',x_noisy);
denoised_wavelet = wdencmp('gbl',x_noisy,'db4',4,thr,sorh,keepapp);
time_wavelet = toc;
denoised_signals{6} = denoised_wavelet;
results(6,:) = [calculate_snr(x_clean, denoised_wavelet), calculate_mse(x_clean, denoised_wavelet), ...
                calculate_rmse(x_clean, denoised_wavelet), calculate_corr(x_clean, denoised_wavelet), ...
                calculate_envelope_smoothness(denoised_wavelet), time_wavelet];

% 3.7 RSBAVMD算法

offset = 0.01;
cutpfreq = 0.0015;
beta = 1e-10;
tol = 1e-7;
ce = 0.3;
tic;
[compset, rsbavmd_iter, comp_freqs] = RSBAVMD(x_noisy, fft(x_noisy), [], fs, tol, ce, offset, cutpfreq);
denoised_rsbavmd = sum(compset(1:end-1,:), 1)'; % 转置为列向量
time_rsbavmd = toc;
denoised_signals{7} = denoised_rsbavmd;
results(7,:) = [calculate_snr(x_clean, denoised_rsbavmd), calculate_mse(x_clean, denoised_rsbavmd), ...
                calculate_rmse(x_clean, denoised_rsbavmd), calculate_corr(x_clean, denoised_rsbavmd), ...
                calculate_envelope_smoothness(denoised_rsbavmd), time_rsbavmd];
iter_info.rsbavmd_iter = rsbavmd_iter;
iter_info.rsbavmd_time = time_rsbavmd;


%% ====================== 4. 窗长与步长敏感性分析 ======================

% 5.1 窗长对迭代次数和时长的影响
window_lengths = 500:100:2500;
vmd_iter_window = zeros(size(window_lengths));
vmd_time_window = zeros(size(window_lengths));
rsbavmd_iter_window = zeros(size(window_lengths));
rsbavmd_time_window = zeros(size(window_lengths));

max_scale = 10; % MPE最大尺度
m = 3; % 嵌入维度
tau = 1; % 延迟时间
window_mpe = cell(length(window_lengths), 1);
window_corr = cell(length(window_lengths), 1);
window_comps = cell(length(window_lengths), 1);
window_denoised = cell(length(window_lengths), 1);

for i = 1:length(window_lengths)
    win_len = window_lengths(i);
    x_window = x_noisy(1:win_len);
    x_clean_window = x_clean(1:win_len);
    
    % VMD
    tic;
    [~,~,~,iter] = VMD(x_window, best_alpha, 0, best_K, 0, 1, 1e-7);
    vmd_time_window(i) = toc;
    vmd_iter_window(i) = iter;
    
    % RSBAVMD
    tic;
    [compset, iter, ~] = RSBAVMD(x_window, fft(x_window), [], fs, tol, ce, offset, cutpfreq);
    rsbavmd_time_window(i) = toc;
    rsbavmd_iter_window(i) = iter;
    
    num_comps = size(compset, 1);
    comp_mpe = zeros(num_comps, max_scale);
    comp_corr = zeros(num_comps, 1);
    
    for j = 1:num_comps
        comp_mpe(j,:) = multiscale_permutation_entropy(compset(j,:)', max_scale, m, tau);
        comp_corr(j) = calculate_corr(x_clean_window, compset(j,:)');
    end
    
    window_mpe{i} = comp_mpe;
    window_corr{i} = comp_corr;
    window_comps{i} = compset;
    window_denoised{i} = sum(compset(1:end-1,:), 1)';
end


% 5.2 步长对迭代次数和时长的影响
step_sizes = 50:10:200;
vmd_iter_step = zeros(size(step_sizes));
vmd_time_step = zeros(size(step_sizes));
rsbavmd_iter_step = zeros(size(step_sizes));
rsbavmd_time_step = zeros(size(step_sizes));

for i = 1:length(step_sizes)
    step = step_sizes(i);
    total_vmd_time = 0;
    total_rsbavmd_time = 0;
    total_vmd_iter = 0;
    total_rsbavmd_iter = 0;
    
    x_current = x_noisy(1:1000);
    fft_current = fft(x_current);
    
    for j = 1:10
        x_new = x_noisy(400+j*step : 400+(j+1)*step-1);
        
        % VMD
        tic;
        [~,~,~,iter] = VMD(x_current, best_alpha, 0, best_K, 0, 1, 1e-7);
        total_vmd_time = total_vmd_time + toc;
        total_vmd_iter = total_vmd_iter + iter;
        
        % RSBAVMD
        tic;
        [~,iter,~] = RSBAVMD(x_current, fft_current, x_new, fs, tol, ce, offset, cutpfreq);
        total_rsbavmd_time = total_rsbavmd_time + toc;
        total_rsbavmd_iter = total_rsbavmd_iter + iter;
        
        x_current = [x_current(step+1:end); x_new];
        fft_current = recursive_FFT(fft_current, x_current, x_new);
    end
    
    vmd_iter_step(i) = total_vmd_iter / 10;
    vmd_time_step(i) = total_vmd_time / 10;
    rsbavmd_iter_step(i) = total_rsbavmd_iter / 10;
    rsbavmd_time_step(i) = total_rsbavmd_time / 10;
end


%% ====================== 4. 不同信噪比下的性能分析 ======================

snr_range = 0:2:50;
vmd_iter_snr = zeros(size(snr_range));
vmd_time_snr = zeros(size(snr_range));
rsbavmd_iter_snr = zeros(size(snr_range));
rsbavmd_time_snr = zeros(size(snr_range));

% 获取RSBAVMD分量数量
[compset_temp,~,~] = RSBAVMD(x_noisy, fft(x_noisy), [], fs, tol, ce, offset, cutpfreq);
num_comps = size(compset_temp, 1);

vmd_rmse_snr = zeros(length(snr_range), best_K);
rsbavmd_rmse_snr = zeros(length(snr_range), num_comps);

for i = 1:length(snr_range)
    snr = snr_range(i);

    noise = randn(size(x_clean));
    noise = noise / norm(noise) * norm(x_clean) / 10^(snr/20);
    x_noisy_snr = x_clean + noise;

    
    % VMD
    tic;
    [u_vmd_snr,~,omega_vmd_snr,iter] = VMD(x_noisy_snr, best_alpha, 0, best_K, 0, 1, 1e-7);
    vmd_time_snr(i) = toc;
    vmd_iter_snr(i) = iter;
    
    % 计算各分量RMSE
    for j = 1:best_K
        vmd_rmse_snr(i,j) = calculate_rmse(x_clean, u_vmd_snr(j,:)');
    end
    
    % RSBAVMD
    tic;
    [compset_snr,iter,~] = RSBAVMD(x_noisy_snr, fft(x_noisy_snr), [], fs, tol, ce, offset, cutpfreq);
    rsbavmd_time_snr(i) = toc;
    rsbavmd_iter_snr(i) = iter;
    
    % 计算各分量RMSE
    for j = 1:min(size(compset_snr,1), num_comps)
        rsbavmd_rmse_snr(i,j) = calculate_rmse(x_clean, compset_snr(j,:)');
    end
end



%% ====================== 5. 真实涡振信号对比实验 ======================

% 加载真实涡振信号（请根据实际路径修改）
data = load('G:\0 新\0 代码（新）\2 涡振噪声信号mat\data_origin_4.mat'); 
Sig = data.data_origin_4(1001:11000,2);

x_real_noisy = Sig;
N_real = length(x_real_noisy);
t_real = (0:N_real-1)/fs;

% 运行所有对比算法
real_methods = {'EMD', 'EEMD', 'Bandpass Filter', 'Grid-VMD', 'Wavelet', 'RSBAVMD'};
num_real_methods = length(real_methods);
real_results = zeros(num_real_methods, 7); % [SNR, MSE, RMSE, Correlation, Envelope Smoothness, SNR Improvement, Time(s)]
real_denoised = cell(num_real_methods, 1);

% 计算原始含噪信号的指标作为基准
base_snr = calculate_snr(x_real_noisy, x_real_noisy); % 原始信号SNR
base_smooth = calculate_envelope_smoothness(x_real_noisy);

% 7.1 EMD算法
tic;
imf_emd_real = emd(x_real_noisy, 'Interpolation','pchip');
denoised_emd_real = sum(imf_emd_real(:,1:4), 2);
time_emd_real = toc;
real_denoised{1} = denoised_emd_real;
real_results(1,1) = calculate_snr(x_real_noisy, denoised_emd_real);
real_results(1,2) = calculate_mse(x_real_noisy, denoised_emd_real);
real_results(1,3) = calculate_rmse(x_real_noisy, denoised_emd_real);
real_results(1,4) = calculate_corr(x_real_noisy, denoised_emd_real);
real_results(1,5) = calculate_envelope_smoothness(denoised_emd_real);
real_results(1,6) = real_results(1,1) - base_snr; % SNR提高量
real_results(1,7) = time_emd_real;

% 7.2 EEMD算法
tic;
imf_eemd_real = eemd_compatible(x_real_noisy, 0.2, 100);
denoised_eemd_real = sum(imf_eemd_real(:,1:4), 2);
time_eemd_real = toc;
real_denoised{2} = denoised_eemd_real;
real_results(2,1) = calculate_snr(x_real_noisy, denoised_eemd_real);
real_results(2,2) = calculate_mse(x_real_noisy, denoised_eemd_real);
real_results(2,3) = calculate_rmse(x_real_noisy, denoised_eemd_real);
real_results(2,4) = calculate_corr(x_real_noisy, denoised_eemd_real);
real_results(2,5) = calculate_envelope_smoothness(denoised_eemd_real);
real_results(2,6) = real_results(2,1) - base_snr;
real_results(2,7) = time_eemd_real;

% 7.3 带通滤波
tic;
[b,a] = butter(4, [0.5 20]/(fs/2), 'bandpass'); % 0.5-20Hz带通滤波
denoised_bp_real = filtfilt(b,a,x_real_noisy);
time_bp_real = toc;
real_denoised{3} = denoised_bp_real;
real_results(3,1) = calculate_snr(x_real_noisy, denoised_bp_real);
real_results(3,2) = calculate_mse(x_real_noisy, denoised_bp_real);
real_results(3,3) = calculate_rmse(x_real_noisy, denoised_bp_real);
real_results(3,4) = calculate_corr(x_real_noisy, denoised_bp_real);
real_results(3,5) = calculate_envelope_smoothness(denoised_bp_real);
real_results(3,6) = real_results(3,1) - base_snr;
real_results(3,7) = time_bp_real;

% 7.4 网格优化VMD
K_range_real = 2:8;
alpha_range_real = 500:500:5000;
best_snr_real = -inf;
best_K_real = 4;
best_alpha_real = 2000;
for K = K_range_real
    for alpha = alpha_range_real
        [u,~,~,~] = VMD(x_real_noisy, alpha, 0, K, 0, 1, 1e-7);
        denoised = sum(u(1:K-1,:), 1)';
        current_snr = calculate_snr(x_real_noisy, denoised);
        if current_snr > best_snr_real
            best_snr_real = current_snr;
            best_K_real = K;
            best_alpha_real = alpha;
        end
    end
end
tic;
[u_vmd_real,~,~,~] = VMD(x_real_noisy, best_alpha_real, 0, best_K_real, 0, 1, 1e-7);
denoised_vmd_real = sum(u_vmd_real(1:best_K_real-1,:), 1)';
time_vmd_real = toc;
real_denoised{4} = denoised_vmd_real;
real_results(4,1) = calculate_snr(x_real_noisy, denoised_vmd_real);
real_results(4,2) = calculate_mse(x_real_noisy, denoised_vmd_real);
real_results(4,3) = calculate_rmse(x_real_noisy, denoised_vmd_real);
real_results(4,4) = calculate_corr(x_real_noisy, denoised_vmd_real);
real_results(4,5) = calculate_envelope_smoothness(denoised_vmd_real);
real_results(4,6) = real_results(4,1) - base_snr;
real_results(4,7) = time_vmd_real;

% 7.5 小波去噪
tic;
[thr_real, sorh_real, keepapp_real] = ddencmp('den','wv',x_real_noisy);
denoised_wavelet_real = wdencmp('gbl',x_real_noisy,'db4',4,thr_real,sorh_real,keepapp_real);
time_wavelet_real = toc;
real_denoised{5} = denoised_wavelet_real;
real_results(5,1) = calculate_snr(x_real_noisy, denoised_wavelet_real);
real_results(5,2) = calculate_mse(x_real_noisy, denoised_wavelet_real);
real_results(5,3) = calculate_rmse(x_real_noisy, denoised_wavelet_real);
real_results(5,4) = calculate_corr(x_real_noisy, denoised_wavelet_real);
real_results(5,5) = calculate_envelope_smoothness(denoised_wavelet_real);
real_results(5,6) = real_results(5,1) - base_snr;
real_results(5,7) = time_wavelet_real;

% 7.6 RSBAVMD算法
tic;
[compset_real,~,~] = RSBAVMD(x_real_noisy, fft(x_real_noisy), [], fs, tol, ce, offset, cutpfreq);
denoised_rsbavmd_real = sum(compset_real(1:end-1,:), 1)';
time_rsbavmd_real = toc;
real_denoised{6} = denoised_rsbavmd_real;
real_results(6,1) = calculate_snr(x_real_noisy, denoised_rsbavmd_real);
real_results(6,2) = calculate_mse(x_real_noisy, denoised_rsbavmd_real);
real_results(6,3) = calculate_rmse(x_real_noisy, denoised_rsbavmd_real);
real_results(6,4) = calculate_corr(x_real_noisy, denoised_rsbavmd_real);
real_results(6,5) = calculate_envelope_smoothness(denoised_rsbavmd_real);
real_results(6,6) = real_results(6,1) - base_snr;
real_results(6,7) = time_rsbavmd_real;


% 滑动窗口分析
window_size = 3000;
step_size = 50;
num_windows = floor((N_real - window_size) / step_size) + 1;

vmd_freqs_real = zeros(num_windows, best_K_real);
vmd_rmse_real = zeros(num_windows, best_K_real);
vmd_iter_real = zeros(num_windows, 1);

rsbavmd_freqs_real = zeros(num_windows, 4); 
rsbavmd_rmse_real = zeros(num_windows, 4);
rsbavmd_iter_real = zeros(num_windows, 1);

for i = 1:num_windows
    start_idx = (i-1)*step_size + 1;
    end_idx = start_idx + window_size - 1;
    x_window = x_real_noisy(start_idx:end_idx);
    
    % VMD分析
    [u_vmd_real,~,omega_vmd_real,iter] = VMD(x_window, best_alpha_real, 0, best_K_real, 0, 1, 1e-7);
    vmd_iter_real(i) = iter;
    for j = 1:best_K_real
        vmd_freqs_real(i,j) = omega_vmd_real(end,j) * fs;
        vmd_rmse_real(i,j) = calculate_rmse(x_window, u_vmd_real(j,:)');
    end
    
    % RSBAVMD分析
    [compset_real,iter,comp_freqs_real] = RSBAVMD(x_window, fft(x_window), [], fs, tol, ce, offset, cutpfreq);
    rsbavmd_iter_real(i) = iter;
    num_comps_real = min(size(compset_real,1), 4);
    for j = 1:num_comps_real
        rsbavmd_freqs_real(i,j) = comp_freqs_real(j);
        rsbavmd_rmse_real(i,j) = calculate_rmse(x_window, compset_real(j,:)');
    end
end


%% ====================== 内嵌函数定义 ======================
% FFT频谱计算函数
function [f, amp] = fft_spectrum(x, fs)
    N = length(x);
    f = (0:N-1)*fs/N;
    f = f(1:floor(N/2)+1);
    Y = fft(x);
    amp = abs(Y/N);
    amp = amp(1:floor(N/2)+1);
    amp(2:end-1) = 2*amp(2:end-1);
end

% 兼容版EEMD函数
function imf = eemd_compatible(x, noise_std, ensemble_num)
    x = x(:);
    N = length(x);
    imf_sum = [];
    
    for i = 1:ensemble_num
        noise = noise_std * randn(size(x));
        x_noisy = x + noise;
        imf_i = emd(x_noisy, 'Interpolation','pchip', 'MaxNumIMF', 10);
        
        if isempty(imf_sum)
            imf_sum = zeros(N, size(imf_i,2));
        end
        
        if size(imf_i,2) < size(imf_sum,2)
            imf_i = [imf_i, zeros(N, size(imf_sum,2)-size(imf_i,2))];
        elseif size(imf_i,2) > size(imf_sum,2)
            imf_sum = [imf_sum, zeros(N, size(imf_i,2)-size(imf_sum,2))];
        end
        
        imf_sum = imf_sum + imf_i;
    end
    
    imf = imf_sum / ensemble_num;
    energy = sum(imf.^2, 1);
    energy_ratio = energy / max(energy);
    imf(:, energy_ratio < 1e-4) = [];
end

% 递归FFT函数
function sig_FFT_new = recursive_FFT(sig_FFT, sig_origin, sig_new_add)
    N = length(sig_origin);
    k = length(sig_new_add);
    x_old_remove = sig_origin(1:k);
    
    u = (0:N-1)';
    shift_factor = exp(-1j * 2 * pi * u * k / N);
    term1 = sig_FFT .* shift_factor;

    term2 = zeros(N, 1);
    for p = 1:k
        delta = sig_new_add(p) - x_old_remove(p);
        if delta ~= 0
            WN = exp(1j * 2 * pi * (k - p + 1) * (0:N-1)' / N);
            term2 = term2 + delta * WN;
        end
    end
    
    sig_FFT_new = term1 + term2;
end

% RSBAVMD主函数
function [compset, total_iter, comp_freqs, error_history] = RSBAVMD(Sig,Sig_FFT,Sig_new_add,SampFreq,tol,ce,offset,cutpfreq)
    Sig1 = Sig;
    N = length(Sig);
    [~,~,~,~,sortInter] = RSpectrendgene(Sig,Sig_FFT,Sig_new_add,SampFreq,offset,cutpfreq);
    maxnum = size(sortInter,1);
    compset = zeros(maxnum,N);
    comp_freqs = zeros(maxnum,1);
    total_iter = 0;
    error_history = [];
    
    
    for ii = 1:maxnum
        Bw = (sortInter(ii,2) - sortInter(ii,1))/SampFreq; 

        % iniIF1 = ((sortInter(ii,2) + sortInter(ii,1))/2)*ones(1,length(Sig));

        iniIF1 = ((sortInter(ii,2) + sortInter(ii,1))/2);

        alpha0 = 500;
        alpha = alpha0*Bw;
        tau = 0;
        DC = 0;
        % init = 1;
        NIMF=1;

        [Sigtemp, ~, omega, iter, error] = RSVMD(Sig, alpha, tau, NIMF, DC, iniIF1, tol);
        Sigtemp = Sigtemp';

        compset(ii,:) = Sigtemp;
        comp_freqs(ii) = omega(end) * SampFreq;
        total_iter = total_iter + iter;
        error_history = [error_history; error];

        Sig = Sig - Sigtemp; 
        coematr = corrcoef(Sigtemp,Sig1);
        if coematr(2) < ce
            break
        end
    end
    compset = compset(1:ii,:); 
    comp_freqs = comp_freqs(1:ii);

    % alpha0 = 500;
    % tau = 0;
    % DC = 0;
    % init = 1;
    % NIMF=maxnum;
    % 
    % [Sigtemp, ~, omega, iter, error] = RSVMD(Sig, alpha0, tau, NIMF, DC, init, tol);
    % Sigtemp = Sigtemp';
    % 
    % compset = Sigtemp;
    % comp_freqs = omega(end) * SampFreq;
    % total_iter = total_iter + iter;
    % error_history = [error_history; error];
    % 
    % Sig = Sig - Sigtemp; 
    % coematr = corrcoef(Sigtemp,Sig1);
    % % if coematr(2) < ce
    % %     break
    % % end
end

% 递归频谱趋势生成函数
function [Spec,Weigt,Spec_trend,WeSpec_trend,sortInter] = RSpectrendgene(Sig,Sig_FFT,Sig_new_add,SampFreq,offset,cutpfreq)
    if isempty(Sig_new_add)
        Sig_FFT_new = fft(Sig);
    else
        Sig_FFT_new = recursive_FFT(Sig_FFT, Sig, Sig_new_add);
    end
    Spec = (2*abs(Sig_FFT_new)/length(Sig));
    Spec = Spec(1:round(end/2));
    Freqbin = linspace(0,SampFreq/2,length(Spec)); 

    orderamp = round(2*length(Sig)*cutpfreq);
    [Spec_trend,~] = coef_ovefour(Spec,SampFreq,orderamp);
    Spec_trend = Spec_trend' + offset;

    [~,indexp] = findpeaks(-Spec_trend);
    tempindex = [1 indexp length(Freqbin)];
    tempindexfrn = [1 indexp];
    tempindexend = [indexp length(Freqbin)];
    diffindex = tempindexend - tempindexfrn;
    deletindex = find(diffindex <= 0.005*length(Sig)); 

    if (~isempty(deletindex)) && (deletindex(1)==1)
        deletindex(1) = deletindex(1) +1;
    end
    tempindex(deletindex) = [];

    Weigt = zeros(size(Freqbin));
    for kk = 1 : length(tempindex)-1
        Weigt(tempindex(kk):tempindex(kk+1)) = sum(Spec(tempindex(kk):tempindex(kk+1)).^2)/(tempindex(kk+1)-tempindex(kk));
    end

    WeSpec_trend = Weigt.*Spec_trend;

    FreqInterval = zeros(length(tempindex)-1,2);
    Weightmax = zeros(1,length(tempindex)-1);
    for dd = 1 : length(tempindex)-1
        Weightmax(dd) = max(WeSpec_trend(tempindex(dd):tempindex(dd+1)-1));
        FreqInterval(dd,:) = [Freqbin(tempindex(dd)),Freqbin(tempindex(dd+1)-1)];
    end

    sortInter = zeros(size(FreqInterval));
    [~,sorin] = sort(Weightmax,'descend');
    for jj = 1:length(tempindex)-1
        sortInter(jj,:) = FreqInterval(sorin(jj),:);
    end
end

% 过完备傅里叶拟合函数
function [fitf,finte] = coef_ovefour(f,SampFreq,orderamp)
    f = f(:);
    N = length(f);
    dt = [0:N-1]/SampFreq;
    f0 = SampFreq/2/N;
    orderamp = 2*orderamp + 1;
    tmatrix = zeros(N,orderamp);
    tmatrix(:,1) = ones(N,1);
    for j = 2:orderamp
        tmatrix(:,j) = cos(2*pi*f0*(j-1)*dt); 
        if j >(orderamp+1)/2
            tmatrix(:,j) = sin(2*pi*f0*(j-((orderamp+1)/2))*dt);
        end
    end
    tmatrix_inte = zeros(N,orderamp);
    tmatrix_inte(:,1) = dt;
    for j = 2:orderamp
        tmatrix_inte(:,j) = 1/(2*pi*f0*(j-1))*sin(2*pi*f0*(j-1)*dt); 
        if j >(orderamp+1)/2
            tmatrix_inte(:,j) = -1/(2*pi*f0*(j-((orderamp+1)/2)))*cos(2*pi*f0*(j-((orderamp+1)/2))*dt);
        end
    end
    alpha = 0.01;
    Imatrix = eye(size(tmatrix,2));
    coeff =(alpha*Imatrix + tmatrix'*tmatrix)\(tmatrix'*f);
    fitf = tmatrix*coeff;
    finte = tmatrix_inte*coeff;
end

% 改进的RSVMD函数（返回迭代次数和误差历史）
function [u, u_hat, omega, iter_num, error_history] = RSVMD(signal, alpha, tau, K, DC, init, tol)
    if isrow(signal)
        signal = signal';
    end
    save_T = length(signal);
    fs = 1/save_T;
    T = save_T;
    f_mirror(1:T/2) = signal(T/2:-1:1);
    f_mirror(1,T/2+1:3*T/2) = signal;
    f_mirror(1,3*T/2+1:2*T) = signal(T:-1:T/2+1);
    f = f_mirror;
    T = length(f);
    t = (1:T)/T;
    freqs = t-0.5-1/T;
    N = 500;
    Alpha = alpha*ones(1,K);
    f_hat = fftshift((fft(f)));
    f_hat_plus = f_hat;
    f_hat_plus(1:T/2) = 0;
    u_hat_plus = zeros(N, length(freqs), K);
    omega_plus = zeros(N, K);
    switch init
        case 1
            for i = 1:K
                omega_plus(1,i) = (0.5/K)*(i-1);
            end
        case 2
            omega_plus(1,:) = sort(exp(log(fs) + (log(0.5)-log(fs))*rand(1,K)));
        otherwise
            omega_plus(1,:) = 0;
    end
    if DC
        omega_plus(1,1) = 0;
    end
    lambda_hat = zeros(N, length(freqs));
    uDiff = tol+eps;
    n = 1;
    sum_uk = 0;
    error_history = [];

    while ( uDiff > tol &&  n < N )
        k = 1;
        sum_uk = u_hat_plus(n,:,K) + sum_uk - u_hat_plus(n,:,1);
        u_hat_plus(n+1,:,k) = (f_hat_plus - sum_uk - lambda_hat(n,:)/2)./(1+Alpha(1,k)*(freqs - omega_plus(n,k)).^2);
        if ~DC
            omega_plus(n+1,k) = (freqs(T/2+1:T)*(abs(u_hat_plus(n+1, T/2+1:T, k)).^2)')/sum(abs(u_hat_plus(n+1,T/2+1:T,k)).^2);
        end
        for k=2:K
            sum_uk = u_hat_plus(n+1,:,k-1) + sum_uk - u_hat_plus(n,:,k);
            u_hat_plus(n+1,:,k) = (f_hat_plus - sum_uk - lambda_hat(n,:)/2)./(1+Alpha(1,k)*(freqs - omega_plus(n,k)).^2);
            omega_plus(n+1,k) = (freqs(T/2+1:T)*(abs(u_hat_plus(n+1, T/2+1:T, k)).^2)')/sum(abs(u_hat_plus(n+1,T/2+1:T,k)).^2);
        end
        lambda_hat(n+1,:) = lambda_hat(n,:) + tau*(sum(u_hat_plus(n+1,:,:),3) - f_hat_plus);
        n = n+1;
        uDiff = eps;
        for i=1:K
            uDiff = uDiff + 1/T*(u_hat_plus(n,:,i)-u_hat_plus(n-1,:,i))*conj((u_hat_plus(n,:,i)-u_hat_plus(n-1,:,i)))';
        end
        uDiff = abs(uDiff);
        error_history = [error_history; uDiff];
    end

    iter_num = n-1;
    omega = omega_plus(1:N,:);
    u_hat = zeros(T, K);
    u_hat((T/2+1):T,:) = squeeze(u_hat_plus(N,(T/2+1):T,:));
    u_hat((T/2+1):-1:2,:) = squeeze(conj(u_hat_plus(N,(T/2+1):T,:)));
    u_hat(1,:) = conj(u_hat(end,:));
    u = zeros(K,length(t));
    for k = 1:K
        u(k,:)=real(ifft(ifftshift(u_hat(:,k))));
    end
    u = u(:,T/4+1:3*T/4);
    clear u_hat;
    for k = 1:K
        u_hat(:,k)=fftshift(fft(u(k,:)))';
    end
end

% 标准VMD函数（返回迭代次数和误差历史）
function [u, u_hat, omega, iter_num, error_history] = VMD(signal, alpha, tau, K, DC, init, tol)
    if isrow(signal)
        signal = signal';
    end
    save_T = length(signal);
    fs = 1/save_T;
    T = save_T;
    f_mirror(1:T/2) = signal(T/2:-1:1);
    f_mirror(1,T/2+1:3*T/2) = signal;
    f_mirror(1,3*T/2+1:2*T) = signal(T:-1:T/2+1);
    f = f_mirror;
    T = length(f);
    t = (1:T)/T;
    freqs = t-0.5-1/T;
    N = 500;
    Alpha = alpha*ones(1,K);
    f_hat = fftshift((fft(f)));
    f_hat_plus = f_hat;
    f_hat_plus(1:T/2) = 0;
    u_hat_plus = zeros(N, length(freqs), K);
    omega_plus = zeros(N, K);
    switch init
        case 1
            for i = 1:K
                omega_plus(1,i) = (0.5/K)*(i-1);
            end
        case 2
            omega_plus(1,:) = sort(exp(log(fs) + (log(0.5)-log(fs))*rand(1,K)));
        otherwise
            omega_plus(1,:) = 0;
    end
    if DC
        omega_plus(1,1) = 0;
    end
    lambda_hat = zeros(N, length(freqs));
    uDiff = tol+eps;
    n = 1;
    sum_uk = 0;
    error_history = [];

    while ( uDiff > tol &&  n < N )
        k = 1;
        sum_uk = u_hat_plus(n,:,K) + sum_uk - u_hat_plus(n,:,1);
        u_hat_plus(n+1,:,k) = (f_hat_plus - sum_uk - lambda_hat(n,:)/2)./(1+Alpha(1,k)*(freqs - omega_plus(n,k)).^2);
        if ~DC
            omega_plus(n+1,k) = (freqs(T/2+1:T)*(abs(u_hat_plus(n+1, T/2+1:T, k)).^2)')/sum(abs(u_hat_plus(n+1,T/2+1:T,k)).^2);
        end
        for k=2:K
            sum_uk = u_hat_plus(n+1,:,k-1) + sum_uk - u_hat_plus(n,:,k);
            u_hat_plus(n+1,:,k) = (f_hat_plus - sum_uk - lambda_hat(n,:)/2)./(1+Alpha(1,k)*(freqs - omega_plus(n,k)).^2);
            omega_plus(n+1,k) = (freqs(T/2+1:T)*(abs(u_hat_plus(n+1, T/2+1:T, k)).^2)')/sum(abs(u_hat_plus(n+1,T/2+1:T,k)).^2);
        end
        lambda_hat(n+1,:) = lambda_hat(n,:) + tau*(sum(u_hat_plus(n+1,:,:),3) - f_hat_plus);
        n = n+1;
        uDiff = eps;
        for i=1:K
            uDiff = uDiff + 1/T*(u_hat_plus(n,:,i)-u_hat_plus(n-1,:,i))*conj((u_hat_plus(n,:,i)-u_hat_plus(n-1,:,i)))';
        end
        uDiff = abs(uDiff);
        error_history = [error_history; uDiff];
    end

    iter_num = n-1;
    omega = omega_plus(1:N,:);
    u_hat = zeros(T, K);
    u_hat((T/2+1):T,:) = squeeze(u_hat_plus(N,(T/2+1):T,:));
    u_hat((T/2+1):-1:2,:) = squeeze(conj(u_hat_plus(N,(T/2+1):T,:)));
    u_hat(1,:) = conj(u_hat(end,:));
    u = zeros(K,length(t));
    for k = 1:K
        u(k,:)=real(ifft(ifftshift(u_hat(:,k))));
    end
    u = u(:,T/4+1:3*T/4);
    clear u_hat;
    for k = 1:K
        u_hat(:,k)=fftshift(fft(u(k,:)))';
    end
end