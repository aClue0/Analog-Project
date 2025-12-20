clearvars;
[x, fs] = audioread('voice.wav');
x = mean(x, 2);
x = (x.')
X = fftshift(fft(x));
len = length(X);
f = (-len/2:len/2-1) * fs / len;
plot(f, abs(X));
xlabel("Frequency (Hz)");
ylabel("|X(f)|");
pause;
close;

filteredX = X .* [zeros(1, floor(len/2) - 4000) , ones(1,8000) , zeros(1, ceil(len/2) - 4000)];

filteredx = real(ifft(ifftshift(filteredX)));
sound (filteredx, fs);
pause;
close;




fc = 100000;
fs2 = 500000;
filteredx = resample(filteredx,fs2,fs);
LEN = length(filteredx)
m = filteredx / max(abs(filteredx));
m = cumsum(m) / fs2;
t = (0:LEN-1) / fs2;
NBFM = 10 * sin(2 * pi * fc * t + 2 * pi * 0.1 * 4000 *m);
WBFM = 10 * sin(2 * pi * fc * t + 2 * pi * 10 * 4000 *m);
#kf = B * fm so 0.1 * 4000 for the narrow case and 10 * 4000 for the wide case
NBFMF = fftshift(fft(NBFM));
WBFMF = fftshift(fft(WBFM));


F = (-LEN/2:LEN/2-1) * fs2 / LEN;
figure;
#NBFM time domain
subplot(2,2,1);
Ns = round(2e-3 * fs2);   % 2 ms window
plot(t(1:Ns), NBFM(1:Ns));
xlabel('Time (s)');
ylabel('Amplitude');
title('NBFM – Time Domain');
grid on;
#WBFM time domain
subplot(2,2,3);
plot(t(1:Ns), WBFM(1:Ns));
xlabel('Time (s)');
ylabel('Amplitude');
title('WBFM – Time Domain');
grid on;
#NBFM spectrum
subplot(2,2,2);
plot(F, abs(NBFMF));
xlabel('Frequency (Hz)');
ylabel('|S(f)|');
title('NBFM – Frequency Domain');
xlim([fc-50000 fc+50000]);
grid on;
#WBFM spectrum
subplot(2,2,4);
plot(F, abs(WBFMF));
xlabel('Frequency (Hz)');
ylabel('|S(f)|');
title('WBFM – Frequency Domain');
xlim([fc-150000 fc+150000]);
grid on;
pause;
close;
#3. What is the condition we needed to achieve NBFM/WBFM?
# we needed beta << 1 so we could negelict it for the narrow case
# we needed beta > 1 so we could NOT negelict it for the wide case



nbfm = diff(NBFM);
wbfm = diff(WBFM);

ned = abs(hilbert(nbfm));
wed = abs(hilbert(wbfm));
ned = ned - mean(ned);
wed = wed - mean(wed);
ned = ned / max(abs(ned));
wed = wed / max(abs(wed));
recnbfm = resample(ned,fs,fs2);
recwbfm = resample(wed,fs,fs2);
Len = length(recnbfm);
recT = (0:Len-1) / fs;
figure;
#NBFM time domain
subplot(2,1,1);
plot(recT, recnbfm);
xlabel('Time (s)');
ylabel('Amplitude');
title('NBFM – Time Domain');
grid on;
#WBFM time domain
subplot(2,1,2);
plot(recT, recwbfm);
xlabel('Time (s)');
ylabel('Amplitude');
title('WBFM – Time Domain');
grid on;
pause;
close;
