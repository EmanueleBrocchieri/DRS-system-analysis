clc;
clear;
close all;

load x.mat
load xclean.mat

figure(1)
plot(x)

fs = 1/0.001; % frequenza di campionamento
freq = fs/length(x)*(0:length(x));
freq = freq(1:floor(length(freq)/2));

plot(x)
title('Segnale rumoroso')
xlabel('Campioni');
ylabel('Ampiezza')
 
% Grafico il modulo dello spettro del segnale
X_mags = abs(fft(x));
figure(10)
plot(X_mags)
xlabel('DFT Bins')
ylabel('Magnitude')
 
% Grafico la prima metà della DFT (frequenze normalizzate)
num_bins = length(X_mags);
plot([0:1/(num_bins/2 -1):1], X_mags(1:num_bins/2))
% plot(freq, X_mags(1:num_bins/2))
xlabel('Frequenza normalizzata (\pi rads/sample)')
ylabel('Magnitude')
 
% Implemento un filtro del secondo ordine utilizzando butterworth 
[b, a] = butter(2, 0.3, 'low');

% b e a sono i coefficienti della funzione di trasferimento del filtro
% H(Z) = (b(1) + b(2)*z^(-1) + ... + b(n+1)*z^(-n))/(a(1) + a(2)*z^(-1) + ... + a(n+1)*z^(-n))
 
% Grafico la risposta in frequenza (frequenza normalizzata)
H = freqz(b,a, floor(num_bins/2));
hold on
plot([0:1/(num_bins/2 -1):1], abs(H),'r');
 
%filtro il segnale utilizzando i coefficienti a e b ottenuti durante
% la progettazione del filtro butterworth
x_filtered = filter(b,a,x);

%Mostro la risposta in frequenza del filtro
figure(100)
freqz(b,a, floor(num_bins/2));

% grafico il segnale filtrato
figure(2)
plot(x_filtered,'r')
title('Segnale Filtrato - Utilizzando Butterworth del Secondo Ordine')
xlabel('Campioni');
ylabel('Ampiezza')
 
% Reimplemento il filtro usando un ordine maggiore
[b2, a2] = butter(20, 0.3, 'low');
 
%Grafico il modulo dello spettro e lo comparo con quello di ordine minore
H2 = freqz(b2,a2, floor(num_bins/2));
figure(10)
hold on
plot([0:1/(num_bins/2 -1):1], abs(H2),'g');
 
% Filtro il segnale e mostro i risultati
x_filtered2 = filter(b2,a2,x);
figure(3)
plot(x_filtered2,'g')
title('Segnale filtrato - Utilizzando Butterworth di ordine 20')
xlabel('Samples');
ylabel('Amplitude')

% Confronto tra segnale pulito e segnale filtrato
figure(200)
plot(xclean,'r')
hold on
plot(x_filtered2,'g')
title('Comparazione tra segnale originario pulito e filtrato')
xlabel('Samples');
ylabel('Amplitude')
legend('originale', 'filtrato');
hold off
 
% Utilizzo un filtro elimina banda al posto del passabasso
[b_stop, a_stop] = butter(20, [0.5 0.8], 'stop');
 
% Grafico il modulo dello spettro
H_stopband = freqz(b_stop,a_stop, floor(num_bins/2));
figure(10)
hold on
plot([0:1/(num_bins/2 -1):1], abs(H_stopband),'c');
 
% Grafico il segnale filtrato
x_filtered_stop = filter(b_stop,a_stop,x);
figure(4);
plot(x_filtered_stop,'c')
title('Segnale Filtrato - Eliminabanda')
xlabel('Campioni');
ylabel('Ampiezza')
 
%Utilizzo la funzione buttord per avere l'ordine ottimale
% Miei requisiti:
% - passband fino a 0.1 (frequenza normalizzata)
% - stopband da 0.5 (frequenza normalizzata)
% - garantire un'attenuazione di non più di 5 dB nella passband
% - garantire un'attenuazione di almeno 40 dB nella stopband
[N, Wn] = buttord(0.1, 0.5, 5, 40);
 
% Utilizzo N e Wn ottenuti sopra per implementare il filtro nel solito modo
[b3, a3] = butter(N, Wn, 'low');

% freqz(b3,a3, floor(num_bins/2))
 
% Grafico il modulo dello spettro
H3 = freqz(b3,a3, floor(num_bins/2));
figure(10);
hold on
plot([0:1/(num_bins/2 -1):1], abs(H3),'k');
figure(10)
 
% filtro il segnale e mostro l'output del filtro
x_filtered3 = filter(b3,a3,x);
figure(5);
plot(x_filtered3,'k')
title(['Segnale Filtrato - Usando un Butterworth di ordine ' num2str(N) ])
xlabel('Campioni');
ylabel('Ampiezza')

% Comparazione con altre tecniche di implementazione (chebyshev e elliptical)
[b_butter, a_butter] = butter(4, 0.2, 'low');
H_butter = freqz(b_butter, a_butter);
 
[b_cheby, a_cheby] = cheby1(4, 0.5, 0.2, 'low');
H_cheby = freqz(b_cheby, a_cheby);
 
[b_ellip, a_ellip] = ellip(4, 0.5, 40, 0.2, 'low');
H_ellip = freqz(b_ellip, a_ellip);
 
% Grafico ogni filtro per compararli
figure(11)
norm_freq_axis = [0:1/(512 -1):1];
plot(norm_freq_axis, abs(H_butter))
hold on
plot(norm_freq_axis, abs(H_cheby),'r')
plot(norm_freq_axis, abs(H_ellip),'g')
legend('Butterworth', 'Chebyshev', 'Elliptical')
xlabel('Frequenza Normalizzata');
ylabel('Modulo')
 
% Grafico in dB in modo da verificare il soddisfacimento delle specifiche
figure(12);
plot(norm_freq_axis, 20*log10(abs(H_butter)))
hold on
plot(norm_freq_axis, 20*log10(abs(H_cheby)),'r')
plot(norm_freq_axis, 20*log10(abs(H_ellip)),'g')
legend('Butterworth', 'Chebyshev', 'Elliptical')
xlabel('Frequenza Normalizzata ');
ylabel('Modulo (dB)')
