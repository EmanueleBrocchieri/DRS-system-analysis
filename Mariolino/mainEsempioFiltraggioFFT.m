clc;
clear;
close all;

% Creiamo un semplice segnale con 2 frequenze
dt = 0.001; % intervallo di campionamento [s]
t = 0:dt:1;
fclean = sin(2*pi*50*t) + sin(2*pi*120*t); % segnale pulito costituito
                                           % dalla somma di 2 frequenze
                                           
f = fclean + 2.5 * randn(size(t)); % aggiungo del rumore

figure(1)
plot(t, fclean, 'LineWidth', 3);
hold on
plot(t, f);
xlabel('tempo [s]');
ylabel('Amplitude [V]');
ylim([-10 10]);

% Calcoliamo la fft del segnale rumoroso
n = length(t);
fhat = fft (f, n);

PSD = fhat .* conj(fhat)/n; % PSD = abs(fhat).^2
freq = 1/(dt*n)*(0:n); % creo il vettore delle frequenze [Hz]
L = 1:floor((n/2));

figure
plot(freq(L), PSD(L));
xlabel('freq [Hz]');
ylabel('PSD [V^2]');
title('PSD del segnale rumoroso')

% Utilizziamo la fft per filtrare il segnale
indices = (PSD > 100); % trovo tutte le frequenze che hanno una potenza
                       % significativa
                        
PSDclean = PSD.*indices;

fhat = indices.*fhat; % poniamo a 0 tutti i coefficienti della fft con
                       % piccola ampiezza
                       
ffilt = ifft(fhat);

figure
plot(t, ffilt);
legend ('Segnale Filtrato');
xlabel('tempo [s]');
ylabel('Ampiezza [V]');
ylim([-10 10]);

figure
plot(t, fclean, '-.', 'LineWidth',3);
hold on
plot(t, ffilt, '-o', 'LineWidth',1);
legend ('Segnale Pulito', 'Segnale Filtrato');
xlabel('tempo [s]');
ylabel('Ampiezza [V]');
ylim([-10 10]);
title('Confronto tra segnale pulito e segnale filtrato')










