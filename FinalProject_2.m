% Digital Communications Project: Performance of Different Modulation Types

% Clear workspace and close all figures
clear all;
close all;
clc ;

% Simulation parameters
numBits = 1e5; % Number of bits
SNR_range = 0:4:60; % SNR range in dB


% Initialize BER matrix
BER_OOK = zeros(size(SNR_range));
BER_PRK = zeros(size(SNR_range));
BER_FSK = zeros(size(SNR_range));


% Generate random binary data vector
data = randi([0 1], 1, numBits);


% Modulation
% OOK modulation
modulated_OOK = data; % No change in bits
    
% PRK modulation
modulated_PRK = 2*data - 1; % Map 0 to -1, 1 to 1
    
% FSK modulation
modulated_FSK = zeros(1, numBits);
for m = 1:numBits
    if data(m) == 0
        modulated_FSK(m) = 1; %exp(1j*0); % Modulate 0 on carrier 1
    else
        modulated_FSK(m) = 1i; %exp(1j*pi/2); % Modulate 1 on orthogonal carrier
    end
end
    



scatterplot(modulated_OOK);
title("Manually Modulated OOK");
scatterplot(modulated_PRK);
title("Manually Modulated PRK");
scatterplot(modulated_FSK);
title("Manually Modulated FSK");







% Calculate transmitted signal power


%Method 2
E_modulated_OOK_beforenoise = sum(abs(modulated_OOK).^2);
power_modulated_OOK_beforenoise = E_modulated_OOK_beforenoise/length(modulated_OOK);

E_modulated_PRK_beforenoise = sum(abs(modulated_PRK).^2);
power_modulated_PRK_beforenoise = E_modulated_PRK_beforenoise/length(modulated_PRK);

E_modulated_FSK_beforenoise = sum(abs(modulated_FSK).^2);
power_modulated_FSK_beforenoise = E_modulated_FSK_beforenoise/length(modulated_FSK);
disp('Transmitted Signal Power before adding noise: ');
disp(['OOK: ', num2str(power_modulated_OOK_beforenoise),', PRK: ',num2str(power_modulated_PRK_beforenoise),' FSK: ',num2str(power_modulated_FSK_beforenoise)]);




% Loop over SNR values---------> Manually
for idx = 1:length(SNR_range)
    
    
    % Calculate SNR from DB to Normal
    snr = 10^(SNR_range(idx)/10);
    % Calculate noise variance
    noise_power = 1/snr;
    
    % Apply noise
    
    
    noise_OOK = sqrt(noise_power/2) * (randn(1, numBits) + 1j * randn(1, numBits) );
    received_OOK = modulated_OOK + noise_OOK;
    
    
    noise_PRK = sqrt(noise_power/2) * (randn(1, numBits) + 1j * randn(1, numBits) );
    received_PRK = modulated_PRK + noise_PRK;
    
    
    
    
    noise_FSK = sqrt(noise_power/2) * (randn(1, numBits) + 1j * randn(1, numBits) );
    received_FSK = modulated_FSK + noise_FSK;
    
    
    %Power after noise
    E_modulated_OOK_afternoise = sum(abs(received_OOK).^2);
    power_modulated_OOK_afternoise = E_modulated_OOK_afternoise/length(received_OOK);
    
    E_modulated_PRK_afternoise = sum(abs(received_PRK).^2);
    power_modulated_PRK_afternoise = E_modulated_PRK_afternoise/length(received_PRK);
    
    E_modulated_FSK_afternoise = sum(abs(received_FSK).^2);
    power_modulated_FSK_afternoise = E_modulated_FSK_afternoise/length(received_FSK);
    
    disp(['SNR: ', num2str(SNR_range(idx)), ' dB']);
    disp(['Power after adding noise: ', ' OOK: ',num2str(power_modulated_OOK_afternoise), ' PRK: ',num2str(power_modulated_PRK_afternoise),' FSK: ',num2str(power_modulated_FSK_afternoise)]);
    
    
    
    
    % Decision
    detected_OOK = real(received_OOK) >= 0.5;%Approx 0.5 % OOK
    detected_PRK = real(received_PRK) >= 0;%Approx 0 % PRK
    detected_FSK = real(received_FSK) < imag(received_FSK); %in FSK if real > imaginary component then it equal zero else equals one

    
    % Calculate number of errors
    errors_OOK = sum(xor(data, detected_OOK));
    errors_PRK = sum(xor(data, detected_PRK));
    errors_FSK = sum(xor(data, detected_FSK));
    
    
    % Calculate Bit Error Rate (BER)
    BER_OOK(idx) = errors_OOK / numBits;
    BER_PRK(idx) = errors_PRK / numBits;
    %BER_FSK(idx) = errors_FSK / numBits;
    [~,BER_FSK(idx)] = biterr(data,detected_FSK) ;

    
    
    scatterplot(received_OOK);
    scatterplot(received_PRK);
    scatterplot(received_FSK);
    
    
end

% Plot BER curve against SNR
figure;
semilogy(SNR_range, BER_OOK, 'bo-', 'LineWidth', 2);
hold on;
semilogy(SNR_range, BER_PRK, 'rx-', 'LineWidth', 2);
semilogy(SNR_range, BER_FSK, 'g+-', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
legend('OOK', 'PRK', 'FSK'); % Add legend here
title('BER vs SNR for Different Modulation Types');



% Evaluate BER curves using MATLAB built-in functions for PSK and PAM modulation

%modulation using modem builtin functions


modem_modulated_OOK = genqammod(data,[0  1]);
scatterplot(modem_modulated_OOK);
title("Built-in OOK");

N = 2; 
modem_modulated_PRK = pskmod(data, 2); 
scatterplot(modem_modulated_PRK);
title("Built-in PRK");



modem_modulated_FSK = genqammod(data,[1  1i]); 
scatterplot(modem_modulated_FSK);
title("Built-in FSK");

% Initialize BER matrix
BER_OOK_modem = zeros(size(SNR_range));
BER_PRK_modem = zeros(size(SNR_range));
BER_FSK_modem = zeros(size(SNR_range));


% Loop over SNR values ---------> using modem-functions
for idx_2 = 1:length(SNR_range)

    % Calculate SNR from DB to Normal
    snr_2 = 10^(SNR_range(idx_2)/10);
    % Calculate noise variance
    noise_variance_2 = 1/sqrt(snr_2);
    % Calculate noise variance
    noise_power_2 = 1/snr_2;
    
    % Apply noise
    
    
    noise_OOK_modem = sqrt(noise_power_2/2) * (randn(1, numBits) + 1j * randn(1, numBits) );
    modem_received_OOK = modem_modulated_OOK + noise_OOK_modem;
    
    
    noise_PRK_modem = sqrt(noise_power_2/2) * (randn(1, numBits) + 1j * randn(1, numBits) );
    modem_received_PRK = modem_modulated_PRK + noise_PRK_modem;
    
    
    noise_FSK_modem = sqrt(noise_power_2/2) * (randn(1, numBits) + 1j * randn(1, numBits) );
    modem_received_FSK = modem_modulated_FSK + noise_FSK_modem;
    
    
    
    modem_OOK_detected=genqamdemod(modem_received_OOK,[0  1]);
    modem_PRK_detected=pskdemod(modem_received_PRK, 2);
    modem_FSK_detected=genqamdemod(modem_received_FSK,[1  1i]);
    
    % Calculate number of errors
    errors_OOK_modem = sum(xor(data, modem_OOK_detected));
    errors_PRK_modem = sum(xor(data, modem_PRK_detected));
    errors_FSK_modem = sum(xor(data, modem_FSK_detected));
   
    
    
    % Calculate Bit Error Rate (BER)
    BER_OOK_modem(idx_2) = errors_OOK_modem / numBits;
    BER_PRK_modem(idx_2) = errors_PRK_modem / numBits;
    BER_FSK_modem(idx_2) = errors_FSK_modem / numBits;
    
    
    
    scatterplot(modem_received_OOK);
    scatterplot(modem_received_PRK);
    scatterplot(modem_received_FSK);
    
    
end


% Plot BER curve against SNR
figure;
semilogy(SNR_range, BER_OOK_modem, 'bo-', 'LineWidth', 2);
hold on;
semilogy(SNR_range, BER_PRK_modem, 'rx-', 'LineWidth', 2);
semilogy(SNR_range, BER_FSK_modem, 'g+-', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
legend('OOK', 'PRK', 'FSK'); % Add legend here
title('BER vs SNR for Different Modulation Types Using built in functions');






% Evaluate the probability of error of 16QAM modulation
%----> First using built in functions
data_16=randi([0 15], 1, numBits);

QAM16_modulatedSignal = qammod(data_16, 16);
scatterplot(QAM16_modulatedSignal);
title("16 QAM");


BER_QAM_16 = zeros(1, length(SNR_range));

for i = 1:length(SNR_range)
    % Add AWGN
    QAM_16noisySignal = awgn(QAM16_modulatedSignal, SNR_range(i), 'measured');

    
    % Demodulate received signal
    receivedData = qamdemod(QAM_16noisySignal, 16);
    
    scatterplot(QAM_16noisySignal);

    % Calculate bit error rate
    [numErrors,BER_QAM_16(i)] = symerr(data_16, receivedData);
   
end

% Plot BER vs SNR
figure;
semilogy(SNR_range, BER_QAM_16, 'o-');
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs SNR for 16-QAM Modulation');






% End of code
