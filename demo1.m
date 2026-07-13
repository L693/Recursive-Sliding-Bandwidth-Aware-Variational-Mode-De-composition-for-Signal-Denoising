

clear; close all; clc;


%%  1. Simulated Signal Generation

fs   = 50;                        % Sampling frequency (Hz)
T    = 60;                        % Total signal duration (s)
t    = (0 : 1/fs : T - 1/fs)';   % Time vector (column vector)
N    = length(t);                 % Total number of samples

% 1.1  Damped vibration mode
s_d  = 5 * exp(-0.1 * t) .* cos(4 * pi * t + pi/4);

% 1.2  Intermodulation components
s_1  = 1.8 * cos(14 * pi * t);    % 5 Hz
s_2  = 0.8 * cos(30 * pi * t);    % 10 Hz

% Composite clean signal
x_clean = s_d + s_1 + s_2;

% Manual noise addition,
snr_input = 2;    % dB
noise_var = var(x_clean) / (10^(snr_input/10));
x_noisy   = x_clean + sqrt(noise_var) * randn(N, 1);



%%  2. RSBAVMD 

tol_rsb     = 1e-7;
ce_rsb      = 0.3;
offset_rsb  = 0.01;
cutpfreq_rsb = 0.0015;

s_d = s_d(1:end-50);
s_1 = s_1(1:end-50);
s_2 = s_2(1:end-50);
x_clean = x_clean(1:end-50);
x_noisy = x_noisy(1:end-50);
Sig = x_noisy;
t = t(1:end-50);
SigFFT = fft(Sig);

SigNewAdd = x_noisy(end-49:end);



tic_rsb = tic;

[compset_rsb,Spec] = RSBAVMD(Sig, SigFFT, SigNewAdd, fs, ...
                      tol_rsb, ce_rsb, offset_rsb, cutpfreq_rsb);
time_rsb = toc(tic_rsb);
compset_rsb = compset_rsb(:,51:end-50);
x_clean = x_clean(51:end-50);
x_noisy = x_noisy(51:end-50);
s_d = s_d(51:end-50);
s_1 = s_1(51:end-50);
s_2 = s_2(51:end-50);
t = t(1:end-100);
nComp_rsb = size(compset_rsb, 1);
x_den_rsb = sum(compset_rsb, 1)';    



%%  4. VMD

alpha_vmd = 2000;
tau_vmd   = 0;
DC_vmd    = 0;
init_vmd  = 1;
tol_vmd   = 1e-7;
K_vmd     = 5;


try
    tic_v = tic;
    [u_vmd, ~, ~] = VMD(x_noisy', alpha_vmd, tau_vmd, K_vmd, ...
                           DC_vmd, init_vmd, tol_vmd);
    time_vmd = toc(tic_v);


    kept_v = [];
    for k = 1:K_vmd
        if abs(corr(x_clean, u_vmd(k,:)')) > ce_rsb
            kept_v = [kept_v; u_vmd(k,:)];  
        end
    end
    if isempty(kept_v),  kept_v = u_vmd(1,:);  end
    x_den_vmd = sum(kept_v, 1)';


catch ME

    x_den_vmd = x_noisy;
    time_vmd  = 0;
end



%%  5. EMD


try
    tic_e = tic;
    imfs  = emd(x_noisy);
    time_emd = toc(tic_e);

    nIMF = size(imfs, 2);
    kept_imf = [];
    for k = 1:nIMF
        if abs(corr(x_clean, imfs(:,k))) > ce_rsb
            kept_imf = [kept_imf, imfs(:,k)];  %#ok<AGROW>
        end
    end
    if isempty(kept_imf),  kept_imf = imfs(:,1);  end
    x_den_emd = sum(kept_imf, 2);


catch ME

    x_den_emd = x_noisy;
    time_emd  = 0;
end



%%  6. Bandpass Filtering Denoising


try
    tic_b = tic;
    x_den_bp = zeros(N, 1);
    for b = 1:nBands
        fL = sortInter(b,1);
        fH = sortInter(b,2);
        if fH - fL < 0.1, continue; end    
        Wn = [fL, fH] * 2 / fs;
        if Wn(1) <= 0,  Wn(1) = 0.001; end
        if Wn(2) >= 1,  Wn(2) = 0.999; end
        [bb, aa] = butter(4, Wn, 'bandpass');
        x_den_bp = x_den_bp + filtfilt(bb, aa, x_noisy);
    end
    time_bp = toc(tic_b);

catch ME

    x_den_bp = x_noisy;
    time_bp  = 0;
end




%%  7. Metrics Calculation


m_MPE = 4;
t_MPE = 1;
Scales = 5;


% ---- 7.1  RSBAVMD Component Metrics ----
mpeComp  = zeros(nComp_rsb, 1);
snrComp  = zeros(nComp_rsb, 1);
corrComp = zeros(nComp_rsb, 1);
for i = 1:nComp_rsb
    c = compset_rsb(i, :)';
    mpeComp(i)  = MPerm(c, m_MPE, t_MPE, Scales);
    snrComp(i)  = calcSNR(x_clean, c);
    corrComp(i) = corr(x_clean, c);
end

% ---- 7.2  Denoising Method Metrics ----
methodNames  = {'RSBAVMD','VMD','EMD','Bandpass'};
denSignals   = {x_den_rsb, x_den_vmd, x_den_emd, x_den_bp};
methodTimes  = [time_rsb, time_vmd, time_emd, time_bp];


mpeMeth  = zeros(5,1);
snrMeth  = zeros(5,1);
corrMeth = zeros(5,1);

for j = 1:4
    s = denSignals{j};
    mpeMeth(j)  = MPerm(s, m_MPE, t_MPE, Scales);
    snrMeth(j)  = calcSNR(x_clean, s);
    corrMeth(j) = corr(x_clean, s);
end

mpeNoisy  = MPerm(x_noisy, m_MPE, t_MPE, Scales);
snrNoisy  = calcSNR(x_clean, x_noisy);
corrNoisy = corr(x_clean, x_noisy);
mpeClean  = MPerm(x_clean, m_MPE, t_MPE, Scales);





[f_x,f_y]=fft_spectrum(x_noisy,50);
f_x=f_x';

figure
subplot(6,1,1)
plot(t,s_d)
ylabel('Amplitude','FontName','Times New Roman','FontSize',12);
title('Damped mode','FontName','Times New Roman','FontSize',12)
set(gca,'defaultAxesTickLabelInterpreter','latex');
yticklabels(strrep(yticklabels,'-','−'));

subplot(6,1,2)
plot(t,s_1)
ylabel('Amplitude','FontName','Times New Roman','FontSize',12);
title('Component 1','FontName','Times New Roman','FontSize',12)
set(gca,'defaultAxesTickLabelInterpreter','latex');
yticklabels(strrep(yticklabels,'-','−'));

subplot(6,1,3)
plot(t,s_2)
ylabel('Amplitude','FontName','Times New Roman','FontSize',12);
title('Component 2','FontName','Times New Roman','FontSize',12)
set(gca,'defaultAxesTickLabelInterpreter','latex');
yticklabels(strrep(yticklabels,'-','−'));


subplot(6,1,4)
plot(t,x_clean)
ylabel('Amplitude','FontName','Times New Roman','FontSize',12);
title('Simulation signal','FontName','Times New Roman','FontSize',12)
set(gca,'defaultAxesTickLabelInterpreter','latex');
yticklabels(strrep(yticklabels,'-','−'));

subplot(6,1,5)
plot(t,x_noisy)
ylabel('Amplitude','FontName','Times New Roman','FontSize',12);
title('Simulation signal with noise','FontName','Times New Roman','FontSize',12)
set(gca,'defaultAxesTickLabelInterpreter','latex');
yticklabels(strrep(yticklabels,'-','−'));

subplot(6,1,6)
plot(f_x,f_y)
xlim([0 25])
ylabel('Amplitude','FontName','Times New Roman','FontSize',12);
title('Spectrum of signal with noise','FontName','Times New Roman','FontSize',12)
set(gca,'defaultAxesTickLabelInterpreter','latex');
yticklabels(strrep(yticklabels,'-','−'));










figure

imf_num = size(compset_rsb,1);
column_plot = 2;
Specpset = [];
fs = 50;
for z=1:imf_num
    
    subplot(imf_num+1,column_plot,2*z-1)
    plot(t,compset_rsb(z,:))
    xlabel('Time(s)','FontName','Times New Roman','FontSize',12);
    ylabel('Amplitude','FontName','Times New Roman','FontSize',12);
    title(['Mode ' num2str(z)],'FontName','Times New Roman','FontSize',12)
    set(gca,'defaultAxesTickLabelInterpreter','latex');
    yticklabels(strrep(yticklabels,'-','−'));

    subplot(imf_num+1,column_plot,2*z)
    [f_x,f_y]=fft_spectrum(compset_rsb(z,:),fs);
    f_x=f_x'; 
    plot(f_x,f_y,'LineWIdth',1)
    xlim([0 25])
    xlabel('Frequency(Hz)','FontName','Times New Roman','FontSize',12);
    ylabel('Amplitude','FontName','Times New Roman','FontSize',12);
    title(['Spectrum of mode ' num2str(z) ],'FontName','Times New Roman','FontSize',12)
    set(gca,'defaultAxesTickLabelInterpreter','latex');
    yticklabels(strrep(yticklabels,'-','−'));
    Specpset = [Specpset;f_y];


end


sig = sum(compset_rsb);
subplot(imf_num+1,column_plot,2*z+1)
plot(t,sig)
xlabel('Time(s)','FontName','Times New Roman','FontSize',12);
ylabel('Amplitude','FontName','Times New Roman','FontSize',12);
title('Signal with noise','FontName','Times New Roman','FontSize',12)
set(gca,'defaultAxesTickLabelInterpreter','latex');
yticklabels(strrep(yticklabels,'-','−'));

subplot(imf_num+1,column_plot,2*z+2)
plot(f_x,Specpset,'linewidth',2);
xlim([0 25])
xlabel('Frequency(Hz)','FontName','Times New Roman','FontSize',12);
ylabel('Amplitude','FontName','Times New Roman','FontSize',12);
title('Spectrum of signal with noise','FontName','Times New Roman','FontSize',12)
legend('Component 1', 'Component 2', 'Component 3','Component 4','FontSize',6,'FontWeight','bold');  % Add legend
set(gca,'defaultAxesTickLabelInterpreter','latex');
yticklabels(strrep(yticklabels,'-','−'));






figure
subplot(3,1,1)
plot(t,x_clean)
xlabel('Time(s)','FontName','Times New Roman','FontSize',12);
ylabel('Amplitude','FontName','Times New Roman','FontSize',12);
ylim([-10 10])
title('Signal without noise','FontName','Times New Roman','FontSize',12)
set(gca,'defaultAxesTickLabelInterpreter','latex');
yticklabels(strrep(yticklabels,'-','−'));


subplot(3,1,2)
plot(t,x_noisy)
ylim([-10 10])
xlabel('Time(s)','FontName','Times New Roman','FontSize',12);
ylabel('Amplitude','FontName','Times New Roman','FontSize',12);
title('Signal with noise','FontName','Times New Roman','FontSize',12)
set(gca,'defaultAxesTickLabelInterpreter','latex');
yticklabels(strrep(yticklabels,'-','−'));

reconstructed_signal = sum(compset_rsb(1:end-1,:));
subplot(3,1,3)
plot(t,reconstructed_signal)
xlabel('Time(s)','FontName','Times New Roman','FontSize',12);
ylabel('Amplitude','FontName','Times New Roman','FontSize',12);
ylim([-10 10])
title('Denoising signal','FontName','Times New Roman','FontSize',12)
set(gca,'defaultAxesTickLabelInterpreter','latex');
yticklabels(strrep(yticklabels,'-','−'));




%% Functions


%  RSBAVMD: Recursive Sliding-window Band-Adaptive VMD
function [compset,Spec] = RSBAVMD(Sig, SigFFT, SigNewAdd, SampFreq, ...
                           tol, ce, offset, cutpfreq)
% RSBAVMD -- Recursive Sliding-window Band-Adaptive VMD
% Input:
%   Sig       - Time-domain signal [N,1] column vector
%   SigFFT    - Previous frame FFT [N,1]
%   SigNewAdd - New samples [M,1]
% Output:
%   compset   - Extracted component matrix [nModes, N]

SigOrig = Sig;                    
N = length(Sig);

% Spectrum analysis: obtain sorted frequency band intervals
[Spec, Weight, SpecTrend, WeSpecTrend, sortInter] = RSpectrendgene(Sig, SigFFT, SigNewAdd, ...
                                          SampFreq, offset, cutpfreq);
maxnum = size(sortInter, 1);

Bw_total = sum(sortInter(:,2) - sortInter(:,1)) / SampFreq;

Alpha = [];


for ii = 1:maxnum

    fLow  = sortInter(ii, 1);
    fHigh = sortInter(ii, 2);
    Bb = (fHigh - fLow) / SampFreq;     
    beta = 500;   
    alpha = beta * (Bw_total - Bb) / Bw_total;
    Alpha = [Alpha,alpha];

end

tau_v  = 0;
DC_v   = 0;
init_v = 1;
K = maxnum;
bw_c = 0.015;


[Sigtemp, ~, omega] = RSVMD(Sig, SigFFT, SigNewAdd, Alpha, tau_v, K, DC_v, init_v, tol, SampFreq);

MPE = [];
compset = [];
omega_c = [];
Sig_n = zeros(1,length(Sig));

for ii = 1:maxnum

    % coemat = corrcoef(Sigtemp(ii, :), SigOrig);
    m_MPE = 4;
    t_MPE = 1;
    Scales = 3;
    MPE_ii = MPerm(Sigtemp(ii, :), m_MPE, t_MPE, Scales);

    if MPE_ii > 0.6 
        Sig_n = Sig_n + Sigtemp(ii, :);
        continue;
    end

    MPE = [MPE;MPE_ii];
    omega_c = [omega_c;omega(ii)];
    f_low = (omega(ii)*2-bw_c);
    if f_low <= 0
        f_high = (omega(ii)*2+bw_c);
        n = 4; 
        Wn = f_high; 
        [b, a] = butter(n, Wn); 
        Sigtemp(ii, :) = filtfilt(b, a, Sigtemp(ii, :));  
    else
        f_high = (omega(ii)*2+bw_c);
        Wn = [f_low, f_high];
        [bb, aa] = butter(4, Wn, 'bandpass');
        Sigtemp(ii, :) = filtfilt(bb, aa, Sigtemp(ii, :));
    end
    compset = [compset;Sigtemp(ii, :)];
    

end

compset = [compset;Sig_n];

end


% RSpectrendgene: Recursive spectrum trend analysis
function [Spec, Weight, SpecTrend, WeSpecTrend, sortInter] = ...
    RSpectrendgene(Sig, SigFFT, SigNewAdd, SampFreq, offset, cutpfreq)

SigFFTnew = recursive_FFT(Sig, SigFFT, SigNewAdd);
SigFFTnew = (2*abs(SigFFTnew)/length(Sig));
SigFFTnew = SigFFTnew(1:round(end/2));

N    = length(Sig);
Spec = 2 * abs(SigFFTnew) / N;
Spec = Spec(1:round(end/2));
Freqbin = linspace(0, SampFreq/2, length(Spec))';

orderAmp = round(2 * N * cutpfreq);
[SpecTrend, ~] = coef_ovefour(Spec, SampFreq, orderAmp);
SpecTrend = SpecTrend' + offset;

[~, indexP] = findpeaks(-SpecTrend);
tempIdx = [1, indexP, length(Freqbin)];
tempIdxFr = [1, indexP];
tempIdxEn = [indexP, length(Freqbin)];
diffIdx = tempIdxEn - tempIdxFr;
deletIdx = find(diffIdx <= 0.005 * N);
if ~isempty(deletIdx)
    if deletIdx(1) == 1
        deletIdx(1) = deletIdx(1) + 1;
    end
    tempIdx(deletIdx) = [];
end

Weight = zeros(size(Freqbin));
for kk = 1 : length(tempIdx)-1
    seg = tempIdx(kk) : tempIdx(kk+1);
    Weight(seg) = sum(Spec(seg).^2) / (tempIdx(kk+1) - tempIdx(kk));
end

WeSpecTrend = Weight .* SpecTrend;

nIntv = length(tempIdx) - 1;
FreqInterval = zeros(nIntv, 2);
WeightMax    = zeros(1, nIntv);
for dd = 1 : nIntv
    seg = tempIdx(dd) : tempIdx(dd+1)-1;
    WeightMax(dd)    = max(WeSpecTrend(seg));
    FreqInterval(dd,:) = [Freqbin(tempIdx(dd)), Freqbin(tempIdx(dd+1)-1)];
end

sortInter = zeros(size(FreqInterval));
[~, sorin] = sort(WeightMax, 'descend');
for jj = 1 : nIntv
    sortInter(jj, :) = FreqInterval(sorin(jj), :);
end
end


% coef_ovefour: Fourier basis fitting
function [fitF, finte] = coef_ovefour(f, SampFreq, orderAmp)
f = f(:);
Nf = length(f);
dt = (0 : Nf-1)' / SampFreq;

f0   = SampFreq / 2 / Nf;
orderAmp = 2 * orderAmp + 1;

tMat = zeros(Nf, orderAmp);
tMat(:,1) = ones(Nf, 1);
for j = 2 : orderAmp
    tMat(:,j) = cos(2 * pi * f0 * (j-1) * dt);
    if j > (orderAmp+1)/2
        tMat(:,j) = sin(2 * pi * f0 * (j - (orderAmp+1)/2) * dt);
    end
end

tMatInte = zeros(Nf, orderAmp);
tMatInte(:,1) = dt;
for j = 2 : orderAmp
    tMatInte(:,j) = 1/(2*pi*f0*(j-1)) * sin(2*pi*f0*(j-1)*dt);
    if j > (orderAmp+1)/2
        tMatInte(:,j) = -1/(2*pi*f0*(j-(orderAmp+1)/2)) * ...
                         cos(2*pi*f0*(j-(orderAmp+1)/2)*dt);
    end
end

alpha  = 0.01;
Imat   = eye(size(tMat, 2));
coeff  = (alpha * Imat + tMat' * tMat) \ (tMat' * f);
fitF   = tMat * coeff;
finte  = tMatInte * coeff;
end


function [u, u_hat, omega] = RSVMD(Sig, SigFFT, SigNewAdd, Alpha, tau, K, DC, init, tol, SampFreq)

%  Input:
%    Sig       - Previous frame original signal [N,1]
%    SigFFT    - Previous frame original signal FFT [N,1] (unshifted)
%    SigNewAdd - New M samples [M,1], M=0 indicates cold start
%    Alpha     - Bandwidth penalty factor, scalar or [1,K] vector
%    tau       - Dual ascent step size, 0=exact reconstruction
%    K         - Number of modes
%    DC        - Whether to fix first mode at DC (0=no)
%    init      - Center frequency initialization (0=0, 1=uniform, 2=random)
%    tol       - Convergence tolerance
%
%  Output:
%    u     - Time-domain modes [K, N]
%    u_hat - Mode spectra [N, K]
%    omega - Center frequencies [K, 1] (normalized 0~0.5)
%


    Sig       = Sig(:);
    SigNewAdd = SigNewAdd(:);
    N         = length(Sig);               % Frame length
    M         = length(SigNewAdd);          % Sliding step size


if isscalar(Alpha)
    Alpha = Alpha * ones(1, K);
end


if M == 0
    signal = Sig;                       % Cold start
else
    signal = [Sig(M+1:N); SigNewAdd];   % Sliding window update
end


if M == 0
    sigFFT = fft(signal);               
else
    sigFFT = recursive_FFT(Sig, SigFFT, SigNewAdd);
end

Niter  = 500;
freqs  = (0:N-1)' / N;                
halfN  = floor(N/2);                
fs_vmd = 1 / N;

u_hat = zeros(N, K);                    
omega = zeros(1, K);

switch init
    case 1
        for i = 1:K
            omega(i) = (0.5/K) * (i-1);
        end
    case 2
        omega = sort(exp(log(fs_vmd) + (log(0.5)-log(fs_vmd)) * rand(1,K)));
    otherwise
        omega(:) = 0;
end
if DC
    omega(1) = 0;
end

lambda_hat = zeros(N, 1);               

% ADMM 
uDiff  = tol + eps;
n      = 0;
prev_u = zeros(N, K);                   

while (uDiff > tol && n < Niter)
    n = n + 1;
    prev_u = u_hat;

    for ki = 1:K

        other_sum = sum(u_hat, 2) - u_hat(:, ki);
        numerator   = sigFFT - other_sum + lambda_hat / 2;
        denominator = 1 + 2 * Alpha(ki) * (freqs - omega(ki)).^2;
        u_hat(:, ki) = numerator ./ denominator;
        omega(ki) = updateCenterFreq(u_hat(:, ki), freqs, halfN);
    end

    if tau > 0
        mode_sum = sum(u_hat, 2);
        lambda_hat = lambda_hat + tau * (sigFFT - mode_sum);
    end

    diff_sq = sum(abs(u_hat(:) - prev_u(:)).^2);
    uDiff   = sqrt(diff_sq) / N;
end

omega = omega(:);

u = zeros(K, N);
for k = 1:K
    u(k, :) = real(ifft(2*u_hat(:, k)))';
end

u_hat_out = zeros(N, K);
for k = 1:K
    u_hat_out(:, k) = fft(u(k, :)');
end
u_hat = u_hat_out;

end


function sigFFTnew = recursive_FFT(sigOrigin, sigFFT, sigNewAdd)
% Incrementally compute current frame FFT from previous frame FFT and new samples
% Input:
%   sigFFT     - Previous frame FFT [N,1] (unshifted)
%   sigOrigin  - Previous frame time-domain signal [N,1]
%   sigNewAdd  - New M samples [M,1]
% Output:
%   sigFFTnew  - Current frame FFT [N,1] (unshifted)

Nf = length(sigOrigin);
k  = length(sigNewAdd);

xOldRemove = sigOrigin(1:k);
deltaX     = sigNewAdd - xOldRemove;

u = (0 : Nf-1)';
shiftFactor = exp(-1j * 2 * pi * u * k / Nf);
term1 = sigFFT .* shiftFactor;

term2 = zeros(Nf, 1);
for p = 1:k
    delta = sigNewAdd(p) - xOldRemove(p);
    if delta ~= 0
        WN = exp(1j * 2 * pi * (k - p + 1) * u / Nf);
        term2 = term2 + delta * WN;
    end
end

sigFFTnew = term1 + term2;

end

function wk = updateCenterFreq(mode_spec, freqs, halfN)
weighted_sum = 0;
power_sum    = 0;
for i = 0:halfN
    pwr = abs(mode_spec(i+1))^2;
    weighted_sum = weighted_sum + freqs(i+1) * pwr;
    power_sum    = power_sum    + pwr;
end

if power_sum > 1e-30
    wk = weighted_sum / power_sum;
else
    wk = 0;
end
end



% VMD
function [u, u_hat, omega] = VMD(signal, alpha, tau, K, DC, init, tol)

% Period and sampling frequency of input signal
save_T = length(signal);
fs = 1/save_T;

% extend the signal by mirroring
T = save_T;
f_mirror(1:T/2) = signal(T/2:-1:1);
f_mirror(T/2+1:3*T/2) = signal;
f_mirror(3*T/2+1:2*T) = signal(T:-1:T/2+1);
f = f_mirror;

% Time Domain 0 to T (of mirrored signal)
T = length(f);
t = (1:T)/T;

% Spectral Domain discretization
freqs = t-0.5-1/T;

% Maximum number of iterations (if not converged yet, then it won't anyway)
N = 500;

% For future generalizations: individual alpha for each mode
Alpha = alpha*ones(1,K);

% Construct and center f_hat
f_hat = fftshift((fft(f)));
f_hat_plus = f_hat;
f_hat_plus(1:T/2) = 0;

% matrix keeping track of every iterant // could be discarded for mem
u_hat_plus = zeros(N, length(freqs), K);

% Initialization of omega_k
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

% if DC mode imposed, set its omega to 0
if DC
    omega_plus(1,1) = 0;
end

% start with empty dual variables
lambda_hat = zeros(N, length(freqs));

% other inits
uDiff = tol+eps; % update step
n = 1; % loop counter
sum_uk = 0; % accumulator

while ( uDiff > tol &&  n < N ) % not converged and below iterations limit
    
    % update first mode accumulator
    k = 1;
    sum_uk = u_hat_plus(n,:,K) + sum_uk - u_hat_plus(n,:,1);
    
    % update spectrum of first mode through Wiener filter of residuals
    u_hat_plus(n+1,:,k) = (f_hat_plus - sum_uk - lambda_hat(n,:)/2)./(1+Alpha(1,k)*(freqs - omega_plus(n,k)).^2);
    
    % update first omega if not held at 0
    if ~DC
        omega_plus(n+1,k) = (freqs(T/2+1:T)*(abs(u_hat_plus(n+1, T/2+1:T, k)).^2)')/sum(abs(u_hat_plus(n+1,T/2+1:T,k)).^2);
    end
    
    % update of any other mode
    for k=2:K
        
        % accumulator
        sum_uk = u_hat_plus(n+1,:,k-1) + sum_uk - u_hat_plus(n,:,k);
        
        % mode spectrum
        u_hat_plus(n+1,:,k) = (f_hat_plus - sum_uk - lambda_hat(n,:)/2)./(1+Alpha(1,k)*(freqs - omega_plus(n,k)).^2);
        
        % center frequencies
        omega_plus(n+1,k) = (freqs(T/2+1:T)*(abs(u_hat_plus(n+1, T/2+1:T, k)).^2)')/sum(abs(u_hat_plus(n+1,T/2+1:T,k)).^2);
        
    end
    
    % Dual ascent
    lambda_hat(n+1,:) = lambda_hat(n,:) + tau*(sum(u_hat_plus(n+1,:,:),3) - f_hat_plus);
    
    % loop counter
    n = n+1;
    
    % converged yet?
    uDiff = eps;
    for i=1:K
        uDiff = uDiff + 1/T*(u_hat_plus(n,:,i)-u_hat_plus(n-1,:,i))*conj((u_hat_plus(n,:,i)-u_hat_plus(n-1,:,i)))';
    end
    uDiff = abs(uDiff);
    
end

N = min(N,n);
omega = omega_plus(1:N,:);

% Signal reconstruction
u_hat = zeros(T, K);
u_hat((T/2+1):T,:) = squeeze(u_hat_plus(N,(T/2+1):T,:));
u_hat((T/2+1):-1:2,:) = squeeze(conj(u_hat_plus(N,(T/2+1):T,:)));
u_hat(1,:) = conj(u_hat(end,:));

u = zeros(K,length(t));

for k = 1:K
    u(k,:)=real(ifft(ifftshift(u_hat(:,k))));
end

% remove mirror part
u = u(:,T/4+1:3*T/4);

% recompute spectrum
clear u_hat;
for k = 1:K
    u_hat(:,k)=fftshift(fft(u(k,:)))';
end

end

% MPE
function MPE = MPerm(X,m,t,Scale)
MPE=[];
for j=1:Scale
    Xs = Multi(X,j);
    PE = pec(Xs,m,t);
    MPE=[MPE PE];
end
MPE = MPE(1,Scale);
MPE = MPE/(log2(factorial(m)));% Normalization
end

function M_Data = Multi(Data,S)
L = length(Data);
J = fix(L/S);
for i=1:J
    M_Data(i) = mean(Data((i-1)*S+1:i*S));
end
end


function [pe hist] = pec(y,m,t)

ly = length(y);
permlist = perms(1:m);
c(1:length(permlist))=0;
    
 for j=1:ly-t*(m-1)
     [a,iv]=sort(y(j:t:j+t*(m-1)));
     for jj=1:length(permlist)
         if (abs(permlist(jj,:)-iv))==0
             c(jj) = c(jj) + 1 ;
         end
     end
 end

hist = c;
 
c=c(find(c~=0));
p = c/sum(c);
pe = -sum(p .* log(p));

end


% calcSNR
function snrVal = calcSNR(clean, estimate)
% SNR = 10 * log10( var(clean) / var(clean - estimate) )
clean    = clean(:);
estimate = estimate(:);
errorSig = clean - estimate;
snrVal   = 10 * log10(var(clean) / var(errorSig));
end


function [f_x,f_y]=fft_spectrum(x,fs)
n=length(x);
hdata=hilbert(x);
f=fft(hdata);
f_x = (0:fs/(n -1):fs)';
f_y=abs(f)/n;
end


