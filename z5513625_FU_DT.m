% === run_taskB_level2.m  (主脚本) ===
clear; clc; close all;

EbN0_dB = 0:2:10;          % 仿真点
bitsPerTrial = 1e5;        % 每个 Eb/N0 生成多少比特
[ber_awgn,  ber_lap] = deal(zeros(size(EbN0_dB)));

for idx = 1:length(EbN0_dB)
    EbN0 = 10^(EbN0_dB(idx)/10);
    N0   = 1/EbN0;                        % 假设 Eb = 1
    % 1 生成比特 & BPSK
    bits  = randi([0 1], bitsPerTrial, 1);
    txSym = 2*bits - 1;                  % 0→-1, 1→+1
    % 2 通道 1：AWGN
    rx1   = txSym + sqrt(N0/2)*randn(bitsPerTrial,1);
    ber_awgn(idx) = mean(sign(rx1) ~= txSym);
    % 3 通道 2：Laplacian
    noiseLap = laprnd(bitsPerTrial,1,0,sqrt(N0/2));  % 自写函数
    rx2   = txSym + noiseLap;
    ber_lap(idx) = mean(sign(rx2) ~= txSym);
end

% After computing ber_awgn(idx) & ber_lap(idx):
ber_awgn_theo = qfunc(sqrt(2*10.^(EbN0_dB/10)));
ber_lap_theo = 0.5*exp(-sqrt(2*10.^(EbN0_dB/10)));

% 4 绘图
semilogy(EbN0_dB, ber_awgn, 'o-', ...
         EbN0_dB, ber_awgn_theo,'--', ...
         EbN0_dB, ber_lap, 's-', ...
         EbN0_dB, ber_lap_theo,'--');
legend('AWGN Sim','AWGN Theory','Laplacian Sim','Lap Theory','Location','southwest');
xlabel('E_b/N_0 (dB)'); ylabel('BER'); grid on;
title('BPSK BER under Different Noise Models (Level 2)');

% === laprnd.m  (生成拉普拉斯分布随机数) ===
function r = laprnd(m,n,mu,b)
% 使用逆变换：U = rand - 0.5 ， X = mu - b * sign(U) .* log(1-2|U|)
U = rand(m,n) - 0.5;
r = mu - b * sign(U) .* log(1 - 2*abs(U));
end

% 检验 Laplacian 分布
r = laprnd(1e5,1,0,1/sqrt(2));   % 单位方差
figure; histogram(r,100,'Normalization','pdf'); hold on;
x = -5:0.01:5; plot(x,1/sqrt(2)*exp(-sqrt(2)*abs(x)),'LineWidth',2);
legend('Empirical','Theoretical'); title('Laplacian PDF');