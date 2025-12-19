% Read the audio file and plot spectrum
[x, Fs] = audioread('testAurora.mp3');
x = x(:,1); % Use mono channel
N = length(x);
t = (0:N-1)/Fs;

% Frequency axis for original signal
f = (-N/2:N/2-1)*(Fs/N);

% Compute spectrum
X = fftshift(fft(x));

% Plot original spectrum
figure;
plot(f, abs(X));
xlabel('Frequency (Hz)');
ylabel('|X(f)|');
title('Original Audio Spectrum');
grid on;

pause; 
close;

% Band-limited filtering (Ideal Low Pass Filter at 4 kHz)
BW = 4000;  % 4 KHz bandwidth
H = abs(f) <= BW;  % Filter mask in frequency domain

Xf = X .* H(:);  % Apply ideal low pass filter in frequency domain

x_filtered = real(ifft(ifftshift(Xf))); % Inverse FFT to obtain time domain filtered signal

% Plot filtered spectrum
figure;
set(gcf, 'Position', [175, 100, 1200, 600]); % wider figure
subplot(1,2,1);
plot(f, abs(Xf));
xlabel('Frequency (Hz)');
ylabel('|X_{filtered}(f)|');
title('Filtered Spectrum (BW = 4 kHz)');
grid on;

% Plot filtered time domain signal
subplot(1,2,2);
plot(t, x_filtered);
xlabel('Time (s)');
ylabel('Amplitude');
title('Filtered Audio Signal in Time Domain');
grid on;

pause;
sound(x_filtered, Fs);  % Play filtered audio
pause;

% Resampling for modulation
Fc = 100e3;  % Carrier frequency
Fs_new = 5 * Fc;  % New sampling frequency (must be >= 2*Fc)

x_up = resample(x_filtered, Fs_new, Fs);
t_up = (0:length(x_up)-1)/Fs_new;

% Generate carrier signals
carrier = cos(2*pi*Fc*t_up(:));

% Generate DSB-SC signal
dsb_sc = x_up .* carrier;

% Generate DSB-TC signal with DC bias
A = 2 * max(abs(x_up));
dsb_tc = (A + x_up) .* carrier;

% Plot spectra of modulated signals
N2 = length(dsb_sc);
f2 = (-N2/2:N2/2-1)*(Fs_new/N2);

DSB_SC_F = fftshift(fft(dsb_sc));
DSB_TC_F = fftshift(fft(dsb_tc));

figure;
set(gcf, 'Position', [175, 100, 1200, 600]); % wider figure
subplot(1,2,1);
plot(f2, abs(DSB_SC_F));
title('DSB-SC Spectrum');
xlabel('Frequency (Hz)');
grid on;

subplot(1,2,2);
plot(f2, abs(DSB_TC_F));
title('DSB-TC Spectrum');
xlabel('Frequency (Hz)');
grid on;

pause;

% Envelope detection using Hilbert transform
env_sc = abs(hilbert(dsb_sc));
env_tc = abs(hilbert(dsb_tc));

% Downsample envelope signals back to original sampling rate for playback
env_sc_ds = resample(env_sc, Fs, Fs_new);
env_tc_ds = resample(env_tc, Fs, Fs_new);

% Play demodulated signals
sound(env_tc_ds, Fs);  % DSB-TC (amplitude envelope)
pause;
sound(env_sc_ds, Fs);  % DSB-SC (envelope should work only if carrier is coherent)
pause;

% --- Coherent detection for DSB-SC ---
SNRs = [0 10 30];  % SNR levels in dB

for i = 1:length(SNRs)
    % Add noise to the DSB-SC signal
    noisy_signal = awgn(dsb_sc, SNRs(i), 'measured');

    % Local carrier for coherent detection
    local_carrier = cos(2*pi*Fc*t_up(:));
    demod = noisy_signal .* local_carrier;

    % Low-pass filtering in frequency domain
    Demod_F = fftshift(fft(demod));
    Demod_F(abs(f2) > BW) = 0;  % Zero out frequencies outside BW
    demod_out = real(ifft(ifftshift(Demod_F)));

    % Downsample to original sampling rate
    demod_ds = resample(demod_out, Fs, Fs_new);

    % Play the demodulated message
    figure;
    plot((0:length(demod_ds)-1)/Fs, demod_ds);
    title(['Recovered Signal at SNR = ', num2str(SNRs(i)), ' dB']);
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;

    sound(demod_ds, Fs);
    pause;
end

% --- Frequency error in carrier ---
Fc_err = 100.1e3;  % Slight frequency error
local_carrier_err = cos(2*pi*Fc_err*t_up);
demod_err = dsb_sc .* local_carrier_err(:);

% Demodulate with frequency error
Demod_F_err = fftshift(fft(demod_err));
Demod_F_err(abs(f2) > BW) = 0;
demod_out_err = real(ifft(ifftshift(Demod_F_err)));

% Downsample and play
demod_ds_err = resample(demod_out_err, Fs, Fs_new);
figure;
plot((0:length(demod_ds_err)-1)/Fs, demod_ds_err);
title('Demodulated with Frequency Error');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
sound(demod_ds_err, Fs);
pause;

% --- Phase error in carrier ---
theta = deg2rad(20);  % Phase error of 20 degrees
local_carrier_phase = cos(2*pi*Fc*t_up + theta);
demod_phase = dsb_sc .* local_carrier_phase(:);

% Demodulate with phase error
Demod_F_phase = fftshift(fft(demod_phase));
Demod_F_phase(abs(f2) > BW) = 0;
demod_out_phase = real(ifft(ifftshift(Demod_F_phase)));

% Downsample and play
demod_ds_phase = resample(demod_out_phase, Fs, Fs_new);
figure;
plot((0:length(demod_ds_phase)-1)/Fs, demod_ds_phase);
title('Demodulated with Phase Error');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
sound(demod_ds_phase, Fs);
pause;

close all; % closes all the plots and figures after completing