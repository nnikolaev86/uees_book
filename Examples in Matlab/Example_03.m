%% ПРИМЕР 3: Устойчивост при малки смущения на модел от 2-ри ред
clear all; clc; close all;
fprintf('ПРИМЕР 3: Устойчивост при малки смущения на модел от 2-ри ред\n')

%% Входни данни

% Режимни параметри на генератора и системата
U1     = 1.02;
P1     = 0.95;

U2     = 1.0;
theta2 = 0; % ъгълът на системата е нула затова ще се пропуска надолу в уравненията

% Параметри на генератора за модел от 2-ри ред
Sn = 500; % MVA
H = 3.0;
Xd_prim  = 0.3;
Ra       = 0.0;

% Собствени и взаимни проводимости без включен генераторен импеданс
Y11  = 0.1882 -2.6990i ;
Y22  = 0.1802 -2.8672i ;
Y12  = 0.1407 -2.5663i ;
Y21  = Y12;

alpha_11 = pi/2 + angle(Y11);
alpha_22 = pi/2 + angle(Y22);
alpha_12 = pi/2 + angle(Y12);
alpha_21 = pi/2 + angle(Y21);

% Собствени и взаимни проводимости с включен генераторен импеданс
% (изчислени са със скрипт Example_03_compute_admitances.m)
Y11g  = 0.0574 -1.4932i ;
Y22g  = 0.0945 -1.7760i ;
Y12g  = 0.0335 -1.4191i ;
Y21g  = Y12g ;

alpha_11g = pi/2 + angle(Y11g);
alpha_22g = pi/2 + angle(Y22g);
alpha_12g = pi/2 + angle(Y12g);
alpha_21g = pi/2 + angle(Y21g);

% Коеф. за преобразуване от RAD в DEG
toDeg = 180/pi; 

%% a) Установен режим на мрежата
fprintf('\na) Установен режим на мрежата\n')

% Ъгъл на генераторното напрежение
theta_10 = asin( (P1 - (U1^2)*abs(Y11)*sin(alpha_11))/(U1*U2*abs(Y12)) ) + alpha_12;
Ug0 = U1*exp(1i*theta_10);
fprintf('Ug0 = (%0.4f /_%0.4f deg) pu\n', abs(Ug0), angle(Ug0)*toDeg)

% Мощности на генератора и системата
P1 = (U1^2)*abs(Y11)*sin(alpha_11) + U1*U2*abs(Y12)*sin(theta_10 - alpha_12);
Q1 = (U1^2)*abs(Y11)*cos(alpha_11) - U1*U2*abs(Y12)*cos(theta_10 - alpha_12);
fprintf('P1      = %0.4f pu\n', P1)
fprintf('Q1      = %0.4f pu\n', Q1)

Ig0 = (P1 - 1i*Q1) / conj(Ug0);
fprintf('Ig0 = %0.4f %+0.4fi pu \n', real(Ig0), imag(Ig0))

%% b) Установен режим на генераторния модел от 2-ри ред
fprintf('\nb) Установен режим на генераторния модел от 2-ри ред\n')

delta0 = angle(Ug0 + (Ra + 1i*Xd_prim)*Ig0);
fprintf('delta0 = %0.4f deg\n', delta0 * toDeg)

Id0 = real(Ig0)*sin(delta0) - imag(Ig0)*cos(delta0);
Iq0 = real(Ig0)*cos(delta0) + imag(Ig0)*sin(delta0);
fprintf('Id0 = %0.4f pu  Iq0 = %0.4f pu \n', Id0, Iq0)

Ud0 = real(Ug0)*sin(delta0) - imag(Ug0)*cos(delta0);
Uq0 = real(Ug0)*cos(delta0) + imag(Ug0)*sin(delta0);
fprintf('Ud0 = %0.4f pu  Uq0 = %0.4f pu \n', Ud0, Uq0)

Eq0_prim = Uq0 + Ra*Iq0 + Xd_prim*Id0;
fprintf('Eq0_prim = %0.4f deg\n', Eq0_prim)



