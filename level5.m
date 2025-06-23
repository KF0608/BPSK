function zID_LastName_DTB_2025_Corrected
    clc; clear; close all;

    % ---------------------------
    % Simulation Parameters
    % ---------------------------
    N = 1e6;                     % Number of Monte Carlo samples
    gamma_dB = 0:1:15;           % SNR range in dB
    gammaLin = 10.^(gamma_dB/10); % Linear SNR values
    C = 1.2;                     % Rate threshold in bps/Hz
    R = 3;                       % Radius for random deployment (m)
    alpha = 2.2;                 % Path-loss exponent
    K = 5;                       % Rician K-factor
    N0 = 1;                      % Noise variance (assumed)
    
    % Threshold for outage calculation
    outage_thresh = 2^C - 1;

    % Arrays to store simulation results
    berAWGN_sim = zeros(size(gammaLin));
    berLap_sim = zeros(size(gammaLin));
    berRayleigh_sim = zeros(size(gammaLin));
    berRician_sim = zeros(size(gammaLin));
    berRayRand_sim = zeros(size(gammaLin));
    berRicRand_sim = zeros(size(gammaLin));

    poutRayleigh_sim = zeros(size(gammaLin));
    poutRician_sim   = zeros(size(gammaLin));
    poutRayRand_sim  = zeros(size(gammaLin));
    poutRicRand_sim  = zeros(size(gammaLin));

    % ---------------------------
    % Start of Main Simulation Loop
    % ---------------------------
    fprintf('Starting Monte Carlo Simulations...\n');
    for iSNR = 1:length(gammaLin)
        g_lin = gammaLin(iSNR);
        fprintf('Simulating SNR = %.1f dB\n', gamma_dB(iSNR));

        % Generate Transmitted Bits
        bits = randi([0 1], N, 1);
        bpskSymbols = 2*bits - 1;  % BPSK: 0->-1, 1->+1

        % --- Channel Models (d=1m or fixed) ---
        
        % 1) AWGN
        noiseAWGN = sqrt(N0/2)*(randn(N,1) + 1j*randn(N,1));
        rxAWGN = sqrt(g_lin)*bpskSymbols + noiseAWGN;
        decAWGN = real(rxAWGN) > 0;
        berAWGN_sim(iSNR) = sum(decAWGN ~= bits)/N;

        % 2) Laplacian Noise
        u = rand(N,1) - 0.5;
        b = sqrt(N0/2); % Scale parameter for Laplacian PDF in doc
        lapNoise = -b * sign(u) .* log(1 - 2*abs(u));
        rxLap = sqrt(g_lin)*bpskSymbols + lapNoise;
        decLap = rxLap > 0;
        berLap_sim(iSNR) = sum(decLap ~= bits)/N;

        % 3) Rayleigh + AWGN
        hRay = (randn(N,1) + 1j*randn(N,1))/sqrt(2);
        rxRay = sqrt(g_lin)*hRay.*bpskSymbols + noiseAWGN;
        rxRayEq = rxRay ./ hRay; 
        decRay = real(rxRayEq) > 0; 
        berRayleigh_sim(iSNR) = sum(decRay ~= bits)/N;
        hRayMagSq = abs(hRay).^2;
        poutRayleigh_sim(iSNR) = mean(g_lin*hRayMagSq < outage_thresh);

        % 4) Rician + AWGN
        hd = exp(1j*2*pi*rand(N,1)); % LoS component with random phase
        hs = (randn(N,1) + 1j*randn(N,1))/sqrt(2); % Scattered component
        hRic = sqrt(K/(K+1))*hd + sqrt(1/(K+1))*hs;
        rxRic = sqrt(g_lin)*hRic.*bpskSymbols + noiseAWGN;
        rxRicEq = rxRic ./ hRic;
        decRic = real(rxRicEq) > 0;
        berRician_sim(iSNR) = sum(decRic ~= bits)/N;
        hRicMagSq = abs(hRic).^2;
        poutRician_sim(iSNR) = mean(g_lin*hRicMagSq < outage_thresh);

        % --- Channel Models (Random Deployment) ---
        dRand = R*sqrt(rand(N,1)); % Random distance d in [0, R]
        path_loss = dRand.^alpha;
        
        % 5) Randomly deployed users - Rayleigh
        hRay2 = (randn(N,1) + 1j*randn(N,1))/sqrt(2);
        effective_snr_ray = g_lin * abs(hRay2).^2 ./ path_loss;
        % For BER, we need channel estimate for equalization
        h_eff_ray = hRay2 ./ (dRand.^(alpha/2));
        rxRayRand = sqrt(g_lin)*h_eff_ray.*bpskSymbols + noiseAWGN;
        rxRayRandEq = rxRayRand ./ h_eff_ray;
        decRayRand  = real(rxRayRandEq) > 0;
        berRayRand_sim(iSNR) = sum(decRayRand ~= bits)/N;
        poutRayRand_sim(iSNR) = mean(effective_snr_ray < outage_thresh);

        % 6) Randomly deployed users - Rician
        hRic2 = sqrt(K/(K+1))*hd + sqrt(1/(K+1))*hs;
        effective_snr_ric = g_lin * abs(hRic2).^2 ./ path_loss;
        h_eff_ric = hRic2 ./ (dRand.^(alpha/2));
        rxRicRand = sqrt(g_lin)*h_eff_ric.*bpskSymbols + noiseAWGN;
        rxRicRandEq = rxRicRand ./ h_eff_ric;
        decRicRand = real(rxRicRandEq) > 0;
        berRicRand_sim(iSNR) = sum(decRicRand ~= bits)/N;
        poutRicRand_sim(iSNR) = mean(effective_snr_ric < outage_thresh);
    end
    fprintf('Simulations finished.\n');

    % ---------------------------
    % ---------------------------
    % Analytical Calculations (Corrected and Stabilized)
    % ---------------------------
    fprintf('Calculating analytical results...\n');
    outage_thresh = 2^C - 1;
    
    % --- Fixed Distance (d=1) ---
    berAWGN_ana = 0.5*erfc(sqrt(gammaLin));
    berLap_ana  = 0.5*exp(-sqrt(2*gammaLin));
    berRayleigh_ana = 0.5*(1 - sqrt(gammaLin./(1+gammaLin)));
    poutRayleigh_ana = 1 - exp(-outage_thresh./gammaLin);
    
    % --- Rician Analytical Calculations (Manual and Stable Implementation) ---
    berRician_ana = zeros(size(gammaLin));
    poutRician_ana = zeros(size(gammaLin));
    
    fprintf('Calculating Rician analytical results (manual integration)...\n');
    for i = 1:length(gammaLin)
        g_lin = gammaLin(i);
        
        % BER for Rician (d=1) using Eq. (13) and (5)
        integrand_ber = @(x) 0.5*erfc(sqrt(x)) .* rician_pdf(x, g_lin, K);
        berRician_ana(i) = integral(integrand_ber, 0, inf);
        
        % Outage for Rician (d=1) using Eq. (17) and STABLE manual Marcum Q-function
        a_marcum = sqrt(2*K);
        b_marcum = sqrt(2*(K+1)*outage_thresh / g_lin);
        poutRician_ana(i) = 1 - MarcumQ1_stable(a_marcum, b_marcum);
    end
    
    % --- Random Deployment (Averaging over distance) ---
    pdf_d = @(d) 2*d/R^2; % PDF of distance d
    
    berRayRand_ana = zeros(size(gammaLin));
    poutRayRand_ana = zeros(size(gammaLin));
    berRicRand_ana = zeros(size(gammaLin));
    poutRicRand_ana = zeros(size(gammaLin));
    
    fprintf('Calculating Random Deployment analytical results...\n');
    for i = 1:length(gammaLin)
        g_lin = gammaLin(i);
        
        % Rayleigh Random Deployment
        ber_integrand_ray = @(d) 0.5*(1 - sqrt((g_lin./d.^alpha)./(1 + g_lin./d.^alpha))) .* pdf_d(d);
        pout_integrand_ray = @(d) (1 - exp(-outage_thresh ./ (g_lin./d.^alpha))) .* pdf_d(d);
        berRayRand_ana(i) = integral(ber_integrand_ray, 0, R);
        poutRayRand_ana(i) = integral(pout_integrand_ray, 0, R);
        
        % Rician Random Deployment (integrating the manual fixed-distance formulas)
        % IMPORTANT: Use 'ArrayValued',true to prevent vectorization error in nested integral
        pout_integrand_ric = @(d) (1 - MarcumQ1_stable(sqrt(2*K), sqrt(2*(K+1)*outage_thresh ./ (g_lin./d.^alpha)))) .* pdf_d(d);
        poutRicRand_ana(i) = integral(pout_integrand_ric, 0, R, 'ArrayValued', true);
        
        ber_integrand_ric = @(d) integral(@(x) 0.5*erfc(sqrt(x)).*rician_pdf(x, g_lin./d.^alpha, K), 0, 50*(1+g_lin./d.^alpha)) .* pdf_d(d);
        berRicRand_ana(i) = integral(ber_integrand_ric, 0, R, 'ArrayValued', true);
    end
    fprintf('Analytical calculations finished.\n');
    % ---------------------------
    % Plotting
    % ---------------------------
    
    % Plot 1: BER for Noise-only channels
    figure; 
    semilogy(gamma_dB, berAWGN_sim, 'bo', 'MarkerFaceColor', 'b'); hold on;
    semilogy(gamma_dB, berAWGN_ana, 'b-', 'LineWidth', 2);
    semilogy(gamma_dB, berLap_sim, 'gs', 'MarkerFaceColor', 'r');
    semilogy(gamma_dB, berLap_ana, 'r-', 'LineWidth', 2);
    grid on; xlabel('SNR \gamma_b (dB)'); ylabel('Bit Error Rate (BER)');
    title('BER Performance in AWGN and Laplacian Noise');
    legend('AWGN Sim','AWGN Ana','Laplacian Sim','Laplacian Ana','Location','southwest');
    ylim([1e-5 1]);

    % Plot 2: BER for Fading channels
    figure;
    semilogy(gamma_dB, berRayleigh_sim, 'bo'); hold on;
    semilogy(gamma_dB, berRayleigh_ana, 'b-', 'LineWidth', 2);
    semilogy(gamma_dB, berRician_sim, 'ro');
    semilogy(gamma_dB, berRician_ana, 'r-', 'LineWidth', 2);
    semilogy(gamma_dB, berRayRand_sim, 'gs');
    semilogy(gamma_dB, berRayRand_ana, 'g--', 'LineWidth', 2);
    semilogy(gamma_dB, berRicRand_sim, 'ms');
    semilogy(gamma_dB, berRicRand_ana, 'm--', 'LineWidth', 2);
    grid on; xlabel('SNR \gamma_b (dB)'); ylabel('Bit Error Rate (BER)');
    title('BER Performance in Fading Channels');
    legend({'Rayleigh Sim','Rayleigh Ana (d=1)', 'Rician Sim','Rician Ana (d=1)', ...
            'Rayleigh Rand Sim','Rayleigh Rand Ana', 'Rician Rand Sim', 'Rician Rand Ana'}, ...
            'Location','southwest');
    ylim([1e-4 1]);

    % Plot 3: Outage Probability for Fading channels
    figure;
    semilogy(gamma_dB, poutRayleigh_sim, 'bo'); hold on;
    semilogy(gamma_dB, poutRayleigh_ana, 'b-', 'LineWidth', 2);
    semilogy(gamma_dB, poutRician_sim, 'ro');
    semilogy(gamma_dB, poutRician_ana, 'r-', 'LineWidth', 2);
    semilogy(gamma_dB, poutRayRand_sim, 'gs');
    semilogy(gamma_dB, poutRayRand_ana, 'g--', 'LineWidth', 2);
    semilogy(gamma_dB, poutRicRand_sim, 'ms');
    semilogy(gamma_dB, poutRicRand_ana, 'm--', 'LineWidth', 2);
    grid on; xlabel('SNR \gamma_b (dB)'); ylabel('Outage Probability (P_{out})');
    title(['Outage Probability for C_{th} = ' num2str(C) ' bps/Hz']);
    legend({'Rayleigh Sim','Rayleigh Ana (d=1)', 'Rician Sim','Rician Ana (d=1)', ...
            'Rayleigh Rand Sim','Rayleigh Rand Ana', 'Rician Rand Sim', 'Rician Rand Ana'}, ...
            'Location','southwest');
    ylim([1e-4 2]);
end


% Helper Functions
% ---------------------------

% Rician PDF (from Eq. 5) - No changes needed, but included for completeness
function pdf = rician_pdf(x, gamma_b, K)
    if gamma_b == 0
        pdf = zeros(size(x));
        return;
    end
    term1 = (K+1)./gamma_b;
    term2 = exp(-K - (K+1).*x./gamma_b);
    term3 = besseli(0, 2*sqrt(K*(K+1).*x./gamma_b));
    pdf = term1 .* term2 .* term3;
    pdf(isnan(pdf) | isinf(pdf)) = 0;
end

% STABLE manual implementation of Marcum Q-function (First Order)
function q = MarcumQ1_stable(a, b)
    % This version is numerically stable by using the scaled Bessel function
    % to avoid Inf*0 issues.
    % Q1(a,b) = integral from b to inf of t*exp(-(t^2+a^2)/2)*I0(a*t) dt
    if b == inf
        q = 0;
        return;
    end
    % Stable integrand: t * besseli(0, a*t, 1) * exp(-(t-a)^2 / 2)
    integrand = @(t) t .* besseli(0, a.*t, 1) .* exp(-(t-a).^2 / 2);
    q = integral(integrand, b, inf);
end