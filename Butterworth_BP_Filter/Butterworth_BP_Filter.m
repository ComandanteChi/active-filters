%% Clear all
clear vars;
close all;
clc;
%% Initialization of variables
ripUp  = 12;                                                               % ripUp...passband upper limit
ripLow = 10;                                                               % ripLow...passband lower limit
Rp = ripUp - ripLow;                                                       % Rp...passband attenuation
Rs = 20;                                                                   % Rs...stopband attenuation
Wp = 2 * pi * [22000 28000];                                               % Wp...pass-frequency
Ws = 2 * pi * [10000 40000];                                               % Ws...stop-frequency
s = tf('s');                                                               % s = tf('s')...allows you to enter the transfer function 
                                                                           %in the next line in symbolic form rather than 
                                                                           %as numerator and denominator vectors.
magnitude = round(10^(ripUp/20));                                          % magnitude...Ac
Ac_1 = -sqrt(9.604e9/(3.08e4*3.85e4));
Ac_2 = Ac_1;
%%  Designing of Butterworth Type I filter
[N, Wc] = buttord(Wp, Ws, Rp, Rs, 's');                                    % N...order of filter, 2*N = 4-> two order sections required
[B, A] = butter(N, Wc, 'bandpass', 's');                                   % B...numerator polynominal H(s)
                                                                           % A...denominator polynominal H(s)
H = magnitude * tf(B, A);                                                  % H...transfer function
P = pole(H);                                                               % P...poles of TF
H1 = Ac_1*(-P(3)-P(4))*s/( (s - P(3)) * (s - P(4)) );                      % H1...TF of 1st order lowpass filter 
H2 = Ac_2*(-P(1)-P(2))*s/( (s - P(1)) * (s - P(2)) );                      % H2...TF of 2nd order lowpass filter
%% H - adaptation of TF to it's prototype
[H_num, H_den] = tfdata(H,'v');                                            % represents TF's numerator and denominator as two vectors

H_num(3) = H_num(3) / H_den(2*N+1);                                        % numerator devided by (2*N+1)-th vector of denominator
H_den(1:1:2*N+1) = H_den(1:1:2*N+1) / H_den(2*N+1);                        % both numerator and denominator are devided by
                                                                           %(2*N-1)th element of vector of denominator
H = tf(H_num, H_den);                                                      % H...transfer function
%%  H1 - adaptation of TF to it's prototype
[H1_num, H1_den] = tfdata(H1,'v');                                         % represents TF's numerator and denominator as two vectors                                          

H1_num(N) = H1_num(N) / H1_den(N+1);                                       % numerator devided by (N+1)-th vector of denominator
H1_den(1:1:N+1) = H1_den(1:1:N+1) / H1_den(N+1);                           % both numerator and denominator are devided by
                                                                           %(N+1)th element of vector of denominator
H1 = tf(H1_num, H1_den);
%%  H2 - adaptation of TF to it's prototype                                                  
[H2_num, H2_den] = tfdata(H2, 'v');                                        % represents TF's numerator and denominator as two vectors

H2_num(N) = H2_num(N) / H2_den(N+1);                                       % numerator devided by (N+1)-th vector of denominator
H2_den(1:1:N+1) = H2_den(1:1:N+1) / H2_den(N+1);                           % both numerator and denominator are devided by
                                                                           %(N+1)th element of vector of denominator
H2 = tf(H2_num, H2_den);
%%  Filter coefficients
d1 = H1_den(1);
c1 = H1_den(2);
d2 = H2_den(1);
c2 = H2_den(2);
%%  Electronic components
C1 = 1e-9;                                                                 % [C1] = F (chosen)
C2 = C1;                                                                   % [C2] = F (chosen)
R1 = -d1 * 2/(c1*Ac_1*C1);                                                 % [R1] = Ohm
R2 = d1/c1 * 2/C1;                                                         % [R2] = Ohm
R3 = 1/C1 * d1/(2*d1/c1 + c1*Ac_1);                                        % [R3] = Ohm
R4 = -d2/c2 * 2/Ac_2*C2;                                                   % [R4] = Ohm
R5 = d2/c2 * 2/C2;                                                         % [R5] = Ohm
R6 = 1/C2 * d2/(2*d2/c2 + c2*Ac_2);                                        % [R6] = Ohm
%%  E24 row
C1_e24 = 1e-9;
C2_e24 = C1_e24;
R1_e24 = 11390;
R2_e24 = 64700;
R3_e24 = 850;
R4_e24 = 9100;
R5_e24 = 51910;
R6_e24 = 680;
%%  transfer functions with the components from the E24
H1_real = tf([0  -R2_e24*R3_e24*C1_e24/(R1_e24+R3_e24)  0],  [R1_e24*R2_e24*R3_e24*C1_e24^2/(R1_e24+R3_e24)   2*R1_e24*R3_e24*C1_e24/(R1_e24+R3_e24)   1]);
H2_real = tf([0  -R5_e24*R6_e24*C2_e24/(R4_e24+R6_e24)  0],  [R4_e24*R5_e24*R6_e24*C2_e24^2/(R4_e24+R6_e24)   2*R4_e24*R6_e24*C2_e24/(R4_e24+R6_e24)   1]);
H_real = H1_real * H2_real;  
%%  Provisional Plot.
subplot(2, 1, 1);
h = bodeplot(H);
setoptions(h, 'FreqUnits', 'kHz', 'PhaseVisible', 'off', 'XLim', {[1, 1000]}, 'YLim', {[-40, 20]});
title('Provisional Plot');
hold on;
h1 = bodeplot(H1);
setoptions(h1, 'FreqUnits', 'kHz', 'PhaseVisible', 'off', 'XLim', {[1, 1000]}, 'YLim', {[-40, 20]});
hold on;
h2 = bodeplot(H2);
setoptions(h2, 'FreqUnits', 'kHz', 'PhaseVisible', 'off', 'XLim', {[1, 1000]}, 'YLim', {[-40, 20]});
hold on;
grid on;
legend('H', 'H1', 'H2', "location", "southwest");

%%  Final Plot.
subplot(2, 1, 2);
h_real = bodeplot(H_real);
setoptions(h_real, 'FreqUnits', 'kHz', 'PhaseVisible', 'off', 'XLim', {[1, 1000]}, 'YLim', {[-40, 20]});
title('Final Plot');
hold on;
h1_real = bodeplot(H1_real);
setoptions(h1_real, 'FreqUnits', 'kHz', 'PhaseVisible', 'off', 'XLim', {[1, 1000]}, 'YLim', {[-40, 20]});
hold on;
h2_real = bodeplot(H2_real);
setoptions(h2_real, 'FreqUnits', 'kHz', 'PhaseVisible', 'off', 'XLim', {[1, 1000]}, 'YLim', {[-40, 20]});
hold on;
grid on;
legend('Hreal', 'H1real', 'H2real', "location", "southwest");