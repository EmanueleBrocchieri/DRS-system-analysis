clc;
clear;
close all;

fs = 1000; % frequenza di campionamento [Hz]
t = 0:(1/fs):1.5-1/fs;

% frequenze del mio segnale [Hz]
f1 = 20; 
f2 = 40;
f3 = 60;

% segnale in input
x = 3*cos(2*pi*f1*t + 0.2) + 1*cos(2*pi*f2*t - 0.3) + ...
    2*cos(2*pi*f3*t + 2.4);

figure(1)
plot(t,x);
xlabel('t [s]');
ylabel('Amplitude [V]');


% calcolo la fft del segnale
X = fft(x);

% z = a + ib
% |z| = sqrt(a^2 + b^2)
% L(z) = tan-1(b/a);

X_mag = abs(X);

figure(2)
plot(X_mag)
xlabel('bins');
ylabel('Amplitude [V]');

% passo dai "bins" alle frequenze
n = length(X);
freq = fs/n*(0:n);

L = 1:floor(n/2);

figure(2)
plot(freq(L), X_mag(L))
xlabel('freq [Hz]');
ylabel('Amplitude [V]');

% Se voglio ottenere il modulo assoluto devo dividere
% per la metà della lunghezza del vettore fft

X_mag = X_mag / length(L);

figure(3)
plot(freq(L), X_mag(L))
xlabel('freq [Hz]');
ylabel('Amplitude [V]');








