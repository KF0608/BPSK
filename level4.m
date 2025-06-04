%% run_taskB_level4.m  ---------------------------------------------------
%  ELEC9123 Design Task B  – Level 4 (Distinction)  
%  功能：在 Level 3 基础上加入随机部署用户的 Rayleigh + AWGN 通道
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

% 随机部署参数
R = 3;                               % 圆形区域半径 (m)
alpha = 2.2;                         % 路径损耗指数
gamma_b = mean(EbN0_lin);            % 平均SNR (用于理论计算)

% 预分配结果向量
ber_awgn = zeros(size(EbN0_dB));
ber_lap  = zeros(size(EbN0_dB));
ber_ray  = zeros(size(EbN0_dB));
ber_ray_random = zeros(size(EbN0_dB));   % 新增：随机部署
outage_ray = zeros(size(EbN0_dB));
outage_ray_random = zeros(size(EbN0_dB)); % 新增：随机部署

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

    % 1.4 Rayleigh 衰落 + AWGN (单位距离)
    h = sqrt(0.5)*(randn(bitsTrial,1)+1i*randn(bitsTrial,1));   % CN(0,1)
    noiseAWGN = sqrt(N0/2)*(randn(bitsTrial,1)+1i*randn(bitsTrial,1));
    rx_ray = h.*txSym + noiseAWGN;
    
    % --- 相位均衡（coherent detection）
    rx_ray_eq = real(conj(h).*rx_ray);
    ber_ray(k) = mean(sign(rx_ray_eq) ~= txSym);
    
    % 1.5 Outage Probability (单位距离)
    gamma_inst = (abs(h).^2) * EbN0;
    outage_ray(k) = mean(gamma_inst < th_SNRlin);
    
    % 1.6 随机部署用户 + Rayleigh 衰落 + AWGN  【新增】
    % 生成随机距离：d ~ f_d(x) = 2x/R^2, 0≤x≤R
    u_dist = rand(bitsTrial,1);                    % 均匀分布 [0,1]
    d = R * sqrt(u_dist);                         % 距离变换
    
    % Rayleigh 衰落系数 (相同)
    h_random = sqrt(0.5)*(randn(bitsTrial,1)+1i*randn(bitsTrial,1));
    
    % 总信道增益：包含衰落和路径损耗
    h_total = h_random ./ (d.^(alpha/2));         % a_r = |h_s|^2 / d^α
    
    % 接收信号
    noiseAWGN_random = sqrt(N0/2)*(randn(bitsTrial,1)+1i*randn(bitsTrial,1));
    rx_ray_random = h_total.*txSym + noiseAWGN_random;
    
    % 相位均衡
    rx_ray_random_eq = real(conj(h_total).*rx_ray_random);
    ber_ray_random(k) = mean(sign(rx_ray_random_eq) ~= txSym);
    
    % Outage：瞬时SNR = |h_total|^2 * Eb/N0  
    gamma_inst_random = (abs(h_total).^2) * EbN0;
    outage_ray_random(k) = mean(gamma_inst_random < th_SNRlin);
end

%% 2. 理论曲线
ber_awgn_theo = qfunc( sqrt(2*EbN0_lin) );
ber_lap_theo  = 0.5*exp( -sqrt(2*EbN0_lin) );

% Rayleigh 理论 BER (单位距离)
ber_ray_theo  = 0.5*( 1 - sqrt( EbN0_lin ./ (1+EbN0_lin) ) );

% Rayleigh Outage (单位距离)
outage_theo   = 1 - exp( - th_SNRlin ./ EbN0_lin );

% 随机部署理论值 (使用近似公式，具体可查阅文档中的积分表达式)
% 这里使用简化理论值作为参考
ber_ray_random_theo = 0.5*( 1 - sqrt( EbN0_lin ./ (2+EbN0_lin) ) ); % 近似
outage_random_theo  = 1 - exp( - 2*th_SNRlin ./ EbN0_lin );         % 近似

%% 3. 绘图

% 3.1 BER 曲线 (所有通道)
figure(1);
semilogy(EbN0_dB, ber_awgn,'bo-', ...
         EbN0_dB, ber_awgn_theo,'b-', ...
         EbN0_dB, ber_lap, 'rs-', ...
         EbN0_dB, ber_lap_theo,'r--', ...
         EbN0_dB, ber_ray, 'kd-', ...
         EbN0_dB, ber_ray_theo,'k--', ...
         EbN0_dB, ber_ray_random,'go-', ...
         EbN0_dB, ber_ray_random_theo,'g--','LineWidth',1.4);
grid on; xlabel('E_b/N_0 (dB)'); ylabel('BER');
title('BPSK BER: AWGN, Laplacian, Rayleigh, and Random Deployment (Level 4)');
legend({'AWGN Sim','AWGN Th.','Lap Sim','Lap Th.','Ray Sim','Ray Th.',...
        'Ray+Random Sim','Ray+Random Th.'},'Location','southwest');

% 3.2 BER 对比 (仅衰落通道)
figure(2);  
semilogy(EbN0_dB, ber_ray, 'kd-', ...
         EbN0_dB, ber_ray_theo,'k--', ...
         EbN0_dB, ber_ray_random,'go-', ...
         EbN0_dB, ber_ray_random_theo,'g--','LineWidth',1.4);
grid on; xlabel('E_b/N_0 (dB)'); ylabel('BER');
title('BER Comparison: Fixed vs Random User Deployment');
legend({'Rayleigh (d=1m) Sim','Rayleigh (d=1m) Th.',...
        'Rayleigh + Random Sim','Rayleigh + Random Th.'},'Location','southwest');

% 3.3 Outage Probability 对比
figure(3);
semilogy(EbN0_dB, outage_ray,'ms-', ...
         EbN0_dB, outage_theo,'m--', ...
         EbN0_dB, outage_ray_random,'co-', ...
         EbN0_dB, outage_random_theo,'c--','LineWidth',1.4);
grid on; xlabel('E_b/N_0 (dB)');
ylabel(['Outage Probability ( \gamma_{th} = ' num2str(th_SNRdB) ' dB )']);
title('Outage Probability: Fixed vs Random User Deployment');
legend({'Rayleigh (d=1m) Sim','Rayleigh (d=1m) Th.',...
        'Rayleigh + Random Sim','Rayleigh + Random Th.'},'Location','northeast');

%% 4. 结果显示
fprintf('\n=== Level 4 Simulation Results ===\n');
fprintf('Parameters: R=%.1fm, α=%.1f, Bits=%d\n', R, alpha, bitsTrial);
fprintf('At Eb/N0=%.0fdB:\n', EbN0_dB(end));
fprintf('  BER (Rayleigh d=1m):     %.2e\n', ber_ray(end));
fprintf('  BER (Rayleigh+Random):   %.2e\n', ber_ray_random(end));  
fprintf('  Outage (Rayleigh d=1m):  %.2e\n', outage_ray(end));
fprintf('  Outage (Rayleigh+Random): %.2e\n', outage_ray_random(end));

disp('=== Level 4 simulation finished ===');

%% 辅助函数
function r = laprnd(m,n,mu,b)
    %LAPRND  Generate i.i.d. Laplacian random numbers
    %   r = laprnd(m,n,mu,b) returns an m-by-n matrix whose entries are
    %   i.i.d. Laplace(mu,b) random variables, where 'b' is the scale
    %   parameter.  Var(X)=2*b^2 .  For unit variance use b = 1/sqrt(2).
    
    U = rand(m,n) - 0.5;                % uniform(-0.5,0.5)
    r = mu - b .* sign(U) .* log(1-2*abs(U));
end
