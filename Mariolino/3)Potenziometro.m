%% Confronto ideale e reale al crescere di Rp/Rm
clc;
clear all;
close all;

V0 = 10.1; % Tensione di alimentazione
L = 260; % Angolo totale di escursione

theta = 0:20:260;

V = V0 * theta/L;

k = 0.1:0.1:1.0; % k = Rp/Rm

Vr = []; % conterra' le future tensioni reali al variare di k

figure
plot(theta, V);

% V = (V0 * theta/L) ./(1 + k(:,ii) * theta/L .* (1 - theta/L)) modello
% reale
hold on

for ii = 1:length(k)
    Vr = [Vr; (V0 * theta/L) ./(1 + k(:,ii) * theta/L .* (1 - theta/L))];
    plot(theta, Vr(ii,:), '-r')
end

hold off
xlabel('angolo [deg]');
ylabel('Vout [V]');

%% OTTIMIZZAZIONE
% DATI 2005
clc;
clear;
close all;

V0 = 10.1; % Tensione di alimentazione
L = 260; % Angolo totale di escursione

theta = [20   40   60   80   100  120  140  160  180  200  220  240  260];
Vout = [0.58 0.8  1.6  1.86 2.56 3.07 3.77 4.05 4.95 5.5  6.1  8.1  9.9]; 
    
% Plot
figure
plot(theta, Vout, 'or')
xlabel('angolo [deg]');
ylabel('Vout [V]');

y = Vout;

% V = V0 * x/L /(1+p*x/L*(1-x/L))
fun = @(p)((V0 * theta/L) ./(1+p*theta/L.*(1-theta/L))- y);

options = optimset('Algorithm','levenberg-marquardt', 'Display','iter', 'TolFun', 1e-20 );
x0 = 0.11;

k = lsqnonlin(fun,x0,[],[],options);


% plotto
Vest = (V0 * theta/L) ./(1 + k * theta/L .* (1-theta/L));

figure
plot(theta, Vest, '-b');
hold on
plot(theta, y, 'ok');
hold off
xlabel('angolo [deg]');
ylabel('Vout [V]');

% calcolo dei residui
r = (V0 * theta/L) ./(1+k*theta/L.*(1-theta/L))- y;

figure
histfit(r);

figure
normplot(r);
