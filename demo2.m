

clear; close all; clc;


%%  1. Signal Generation

fs   = 50;                        % Sampling frequency (Hz)
T    = 200;                        % Total signal duration (s)
t    = (0 : 1/fs : T - 1/fs)';   % Time vector (column vector)
N    = length(t);                 % Total number of samples

data = load('data.mat'); 
Sig = data.Sig;


[f_x,f_y]=fft_spectrum(Sig,50);
f_x=f_x';

figure
subplot(2,1,1)
plot(t,Sig)
xlim([0 150])
xlabel('Time(s)','FontName','Times New Roman','FontSize',16);
ylabel('Amplitude','FontName','Times New Roman','FontSize',16);
title('Signal')
set(gca,'defaultAxesTickLabelInterpreter','latex');
yticklabels(strrep(yticklabels,'-','−'));

subplot(2,1,2)
plot(f_x,f_y)
xlim([-0.2 25])
xlabel('Frequency(Hz)','FontName','Times New Roman','FontSize',16);
ylabel('Amplitude','FontName','Times New Roman','FontSize',16);
title('Spectrum of signal','FontName','Times New Roman','FontSize',16)
set(gca,'defaultAxesTickLabelInterpreter','latex');
yticklabels(strrep(yticklabels,'-','−'));


%%  2. RSBAVMD 

tol_rsb     = 1e-7;
ce_rsb      = 0.3;
offset_rsb  = 0.01;
cutpfreq_rsb = 0.0015;


Sig = Sig(1:end-50);
t = t(1:end-50);
SigFFT = fft(Sig);

SigNewAdd = Sig(end-49:end);



tic_rsb = tic;

[compset_rsb,Spec] = RSBAVMD(Sig, SigFFT, SigNewAdd, fs, ...
                      tol_rsb, ce_rsb, offset_rsb, cutpfreq_rsb);
time_rsb = toc(tic_rsb);
compset_rsb = compset_rsb(:,51:end-50);
Sig = Sig(51:end-50);
t = t(1:end-100);
nComp_rsb = size(compset_rsb, 1);
x_den_rsb = sum(compset_rsb, 1)';    



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
title('Signal','FontName','Times New Roman','FontSize',12)
set(gca,'defaultAxesTickLabelInterpreter','latex');
yticklabels(strrep(yticklabels,'-','−'));

subplot(imf_num+1,column_plot,2*z+2)
plot(f_x,Specpset,'linewidth',2);
xlim([0 25])
xlabel('Frequency(Hz)','FontName','Times New Roman','FontSize',12);
ylabel('Amplitude','FontName','Times New Roman','FontSize',12);
title('Spectrum of signal','FontName','Times New Roman','FontSize',12)
legend('Component 1', 'Component 2', 'Component 3','Component 4','FontSize',6,'FontWeight','bold');  % 添加图例
set(gca,'defaultAxesTickLabelInterpreter','latex');
yticklabels(strrep(yticklabels,'-','−'));


figure
subplot(2,1,1)
plot(t,Sig)
xlabel('Time(s)','FontName','Times New Roman','FontSize',16);
ylabel('Amplitude','FontName','Times New Roman','FontSize',16);
title('Original signal','FontName','Times New Roman','FontSize',16)
set(gca,'defaultAxesTickLabelInterpreter','latex');
yticklabels(strrep(yticklabels,'-','−'));

reconstructed_signal = sum(compset_rsb(1:end-1,:),1);
subplot(2,1,2)
plot(t,reconstructed_signal)
xlabel('Time(s)','FontName','Times New Roman','FontSize',16);
ylabel('Amplitude','FontName','Times New Roman','FontSize',16);
title('Denoising signal','FontName','Times New Roman','FontSize',16)
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
bw_c = 0.005;

[Sigtemp, ~, omega] = RSVMD(Sig, SigFFT, SigNewAdd, Alpha, tau_v, K, DC_v, init_v, tol, SampFreq);

MPE = [];
compset = [];
omega_c = [];
Sig_n = zeros(1,length(Sig));

for ii = 1:maxnum

    coemat = corrcoef(Sigtemp(ii, :), SigOrig);
    m_MPE = 4;
    t_MPE = 1;
    Scales = 3;
    MPE_ii = MPerm(Sigtemp(ii, :), m_MPE, t_MPE, Scales);

    if MPE_ii > 0.6  || max(Sigtemp(ii, :)) > 0.1*max(Sig)
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


