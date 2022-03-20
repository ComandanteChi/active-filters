%% Clear all
clear vars;
close all;
clc;
%% Initialization of variables
ripUp  = 9.5;                                                               % ripUp...passband upper limit
ripLow = 3.5;                                                               % ripLow...passband lower limit
Rp = ripUp - ripLow;                                                       % Rp...passband attenuation
Rs = 45;                                                                   % Rs...stopband attenuation
Wp = 20000 * 2 * pi;                                                       % Wp...pass-frequency
Ws = 50000 * 2 * pi;                                                       % Ws...stop-frequency
s = tf('s');                                                               % s = tf('s')...allows you to enter the transfer function 
                                                                           %in the next line in symbolic form rather than 
                                                                           %as numerator and denominator vectors.
magnitude = 10^(ripUp/20);                                                 % magnitude...A0
%%  Designing of Chebyshev Type I filter
[N, Wp] = cheb1ord(Wp, Ws, Rp, Rs, 's');                                    % N...order of filter and cuttoff frequency (Grenzfrequenz)
[B, A] = cheby1(N, Rp, Wp, 's');                                           % B...numerator polynominal H(s)
                                                                           % A...denominator polynominal H(s)
H = magnitude * tf(B, A);                                                  % H...transfer function
P = pole(H);                                                               % P...poles of TF
H1 = sqrt(magnitude)*(-P(3))*(-P(4))/( (s - P(3)) * (s - P(4)) );          % H1...TF of 2 order lowpass filter 
H2 = sqrt(magnitude)*(-P(1))*(-P(2))/( (s - P(1)) * (s - P(2)) );          % H2...TF of 2nd order lowpass filter
%% H - adaptation of TF to it's prototype
[H_num, H_den] = tfdata(H,'v')                                            % represents TF's numerator and denominator as two vectors

H_num(N+1)     = H_num(N+1)     / H_den(N+1);                              % numerator devided by (N+1)-th vector of denominator
H_den(1:1:N+1) = H_den(1:1:N+1) / H_den(N+1);                              % both numerator and denominator are devided by
                                                                           %(N+1)th element of vector of denominator
H = tf(H_num, H_den);                                                      % H...transfer function
%%  H1 - adaptation of TF to it's prototype
[H1_num, H1_den] = tfdata(H1,'v');                                         % represents TF's numerator and denominator as two vectors                                          

H1_num(3) = H1_num(3) / H1_den(3);                                         % numerator devided by (N+1)-th vector of denominator
H1_den(1:1:3) = H1_den(1:1:3) / H1_den(3);                                 % both numerator and denominator are devided by
                                                                           %(N+1)th element of vector of denominator
H1 = tf(H1_num, H1_den);
%%  H2 - adaptation of TF to it's prototype                                                  
[H2_num, H2_den] = tfdata(H2, 'v');                                        % represents TF's numerator and denominator as two vectors

H2_num(3) = H2_num(3) / H2_den(3);                                         % numerator devided by (N+1)-th vector of denominator
H2_den(1:1:3) = H2_den(1:1:3) / H2_den(3);                                 % both numerator and denominator are devided by
                                                                           %(N+1)th element of vector of denominator
H2 = tf(H2_num, H2_den);
%%  Filter coefficients
a1 = H1_den(1);
a2 = H2_den(2);
b2 = H2_den(1);
%%  Electronic components
C1 = 10e-9;                                                                % [C1] = F (chosen)
C2 = 680e-12;                                                              % [C2] = F (chosen)
%C3 = ( 4 * b2 * (1 + H2_num(3)) * C2 ) / a2^2;                            % [C3] = F
C3 = 68e-9;                                                                % [C3] = F
R2 = a1/C1;                                                                % Ohm
R1 = R2/H1_num(2);                                                         % Ohm
R4 = ( a2 * C3 - sqrt(a2^2 * C3^2 - 4 * C2 * C3 * b2 * (1 + H2_num(3))) ) / ( 2 * C2 * C3);     % Ohm
R3 = R4 / H2_num(3);                                                       % Ohm 
R5 = b2 / (C2 * C3 * R4);                                                  % Ohm
%%  E24 row
C1_e24 = 10e-9;                                                            % [C1] = F (chosen)
C2_e24 = 680e-12;                                                          % [C2] = F (chosen)
C3_e24 = 68e-9;                                                            % [C3] = F
R2_e24 = 4300;                                                             % Ohm
R1_e24 = 2000 + 160;                                                       % Ohm
R4_e24 = 2200 + 390;                                                       % Ohm
R3_e24 = 1300;                                                             % Ohm
R5_e24 = 2400;                                                             % Ohm
%%  transfer functions with the components from the E24
H1_real = tf([0 (-R2_e24/R1_e24)], [(R2_e24*C1_e24) 1]);
H2_real = tf([0 0 (-R4_e24/R3_e24)], [(C2_e24*C3_e24*R4_e24*R5_e24) (C2_e24*(R4_e24+R5_e24+R4_e24*R5_e24/R3_e24)) 1]);
H_real = H1_real*H2_real;
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
legend('H', 'H1', 'H3', "location", "southwest");

%%  Final Plot.
subplot(2, 1, 2);
h = bodeplot(H_real);
setoptions(h, 'FreqUnits', 'kHz', 'PhaseVisible', 'off', 'XLim', {[1, 1000]}, 'YLim', {[-40, 20]});
title('Final Plot');
hold on;
h1 = bodeplot(H1_real);
setoptions(h1, 'FreqUnits', 'kHz', 'PhaseVisible', 'off', 'XLim', {[1, 1000]}, 'YLim', {[-40, 20]});
hold on;
h2 = bodeplot(H2_real);
setoptions(h2, 'FreqUnits', 'kHz', 'PhaseVisible', 'off', 'XLim', {[1, 1000]}, 'YLim', {[-40, 20]});
hold on;
grid on;
legend('H', 'H1', 'H3', "location", "southwest");
