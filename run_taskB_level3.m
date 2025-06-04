%% run_taskB_level3.m  ---------------------------------------------------
%  ELEC9123 Design Task B  – Level 3 (Credit)
%  功能：在 Level 2 基础上再加入 Rayleigh + AWGN 通道，并计算 Outage Probability
%  作者：<Your Name>  zID:<z530xxxx>   日期：2025-MM-DD
%  运行环境：MATLAB R2023b + Communications Toolbox

clear; clc; close all;
rng(2025);                           % 固定随机种子，方便助教复现

%% 0. 仿真参数
EbN0_dB   = 0:2:10;                  % x 轴
EbN0_lin  = 10.^(EbN0_dB/10);        % 线性刻度
bitsTrial = 1e5;                     % 每个 Eb/N0 仿真的比特数
th_SNRdB  = 10;                      % Outage 门限 γ_th = 10 dB
th_SNRlin = 10^(th_SNRdB/10);

% 预分配结果向量
ber_awgn = zeros(size(EbN0_dB));
ber_lap  = zeros(size(EbN0_dB));
ber_ray  = zeros(size(EbN0_dB));
outage_ray = zeros(size(EbN0_dB));

%% 1. 主循环
for k = 1:length(EbN0_dB)
    EbN0 = EbN0_lin(k);
    N0   = 1/EbN0;                   % 假设 Eb = 1

    % 1.1 生成比特 & BPSK 映射
    bits  = randi([0 1], bitsTrial, 1);
    txSym = 2*bits - 1;              % {0,1} -> {-1,+1}

    % 1.2 AWGN 通道
    rx_awgn = txSym + sqrt(N0/2)*randn(bitsTrial,1);
    ber_awgn(k) = mean(sign(rx_awgn) ~= txSym);

    % 1.3 Laplacian 噪声通道
    noiseLap = laprnd(bitsTrial,1,0,sqrt(N0/2));
    rx_lap   = txSym + noiseLap;
    ber_lap(k) = mean(sign(rx_lap) ~= txSym);

    % 1.4 Rayleigh 衰落 + AWGN
    h = sqrt(0.5)*(randn(bitsTrial,1)+1i*randn(bitsTrial,1));   % CN(0,1)
    noiseAWGN = sqrt(N0/2)*(randn(bitsTrial,1)+1i*randn(bitsTrial,1));
    rx_ray = h.*txSym + noiseAWGN;

    % --- 相位均衡（coherent detection）
    rx_ray_eq = real(conj(h).*rx_ray);
    ber_ray(k) = mean(sign(rx_ray_eq) ~= txSym);

    % 1.5 Outage Probability  γ = |h|^2 · Eb/N0  < γ_th？
    gamma_inst = (abs(h).^2) * EbN0;
    outage_ray(k) = mean(gamma_inst < th_SNRlin);
end

%% 2. 理论曲线
ber_awgn_theo = qfunc( sqrt(2*EbN0_lin) );
ber_lap_theo  = 0.5*exp( -sqrt(2*EbN0_lin) );
% Rayleigh 理论 BER：Pe = 0.5*(1 - sqrt( Eb/N0 ./ (1+Eb/N0) ))
ber_ray_theo  = 0.5*( 1 - sqrt( EbN0_lin ./ (1+EbN0_lin) ) );
% Rayleigh Outage：P_out = 1 - exp( - γ_th / Eb/N0 )
outage_theo   = 1 - exp( - th_SNRlin ./ EbN0_lin );

%% 3. 绘图

% 3.1 BER 曲线
figure(1);
semilogy(EbN0_dB, ber_awgn,'bo-', ...
         EbN0_dB, ber_awgn_theo,'b-', ...
         EbN0_dB, ber_lap, 'rs-', ...
         EbN0_dB, ber_lap_theo,'r--', ...
         EbN0_dB, ber_ray, 'kd-', ...
         EbN0_dB, ber_ray_theo,'k--','LineWidth',1.4);
grid on; xlabel('E_b/N_0 (dB)'); ylabel('BER');
title('BPSK BER under AWGN, Laplacian and Rayleigh Channels (Level 3)');
legend({'AWGN Sim','AWGN Th.','Lap Sim','Lap Th.','Ray Sim','Ray Th.'},...
       'Location','southwest');

% 3.2 Outage Probability
figure(2);
semilogy(EbN0_dB, outage_ray,'ms-','LineWidth',1.4); hold on;
semilogy(EbN0_dB, outage_theo,'m--','LineWidth',1.4);
grid on; xlabel('E_b/N_0 (dB)');
ylabel(['Outage P  ( \gamma_{th} = ' num2str(th_SNRdB) ' dB )']);
title('Rayleigh Channel Outage Probability');
legend({'Simulation','Theory'},'Location','northeast');

disp('=== Level 3 simulation finished ===');

function r = laprnd(m,n,mu,b)
    %LAPRND  Generate i.i.d. Laplacian random numbers
    %   r = laprnd(m,n,mu,b) returns an m-by-n matrix whose entries are
    %   i.i.d. Laplace(mu,b) random variables, where 'b' is the scale
    %   parameter.  Var(X)=2*b^2 .  For unit variance use b = 1/sqrt(2).
    
    U = rand(m,n) - 0.5;                % uniform(-0.5,0.5)
    r = mu - b .* sign(U) .* log(1-2*abs(U));
    end