%% ПРИМЕР 2: Изчисляване на ъглови характеристики и установен режим на едномашинна ЕЕС
clear all; clc; close all;
fprintf('ПРИМЕР 2: Изчисляване на ъглови характеристики и установен режим на едномашинна ЕЕС\n')

%% Входни данни

%Собствени проводимости:
Y11  = 0.1882 -2.6990i ;
Y22  = 0.1802 -2.8672i ;
%Взаимни проводимости:
Y12  = 0.1407 -2.5663i ;
Y21  = Y12;

% За генераторния възел
U1 = 1.02;
P1 = 0.95;
% За системата
U2 = 1.0;
theta2 = 0; % ъгълът на системата е нула затова ще се пропуска надолу в уравненията

% Коеф. за преобразуване от RAD в DEG
toDeg = 180/pi; 

% Параметри на генератора
Sn = 500; % MVA
H = 3.0;
Td0_prim = 6.5; % s
Td0_sec  = 0.035;
Tq0_prim = 1.4;
Tq0_sec  = 0.04;
Xd       = 1.9; % pu
Xq       = 1.8;
Xd_prim  = 0.3;
Xq_prim  = 0.55;
Xd_sec   = 0.26;
Xq_sec   = Xd_sec;
Xl       = 0.2;
Ra       = 0.0;
K1 = (Xd - Xd_prim) / (Xd_prim - Xl);
K2 = 1 + K1;
K3 = (Xd - Xd_prim) * (1 - (Xd_prim - Xd_sec) / (Xd_prim - Xl));

%% а) Ъглови характеристики
fprintf('\na) Ъглови характеристики\n')

% Допълващи ъгли на проводимостите (всички са в RAD)
alpha_11 = pi/2 + angle(Y11);
alpha_22 = pi/2 + angle(Y22);
alpha_12 = pi/2 + angle(Y12);
alpha_21 = pi/2 + angle(Y21);

fprintf('alpha_11 = %0.4f deg\n', alpha_11 * toDeg)
fprintf('alpha_22 = %0.4f deg\n', alpha_22 * toDeg)
fprintf('alpha_12 = %0.4f deg\n', alpha_12 * toDeg)
fprintf('alpha_21 = %0.4f deg\n\n', alpha_21 * toDeg)

% Диапазон на изменение на ъгъл theta_1
theta_1 = linspace(0, pi, 100); % генериране на 100 точки между 0 и pi/2 RAD (0 до 90 DEG)

% Ъглови х-ки на генератора
P1_theta = (U1^2)*abs(Y11)*sin(alpha_11) + U1*U2*abs(Y12)*sin(theta_1 - alpha_12);
Q1_theta = (U1^2)*abs(Y11)*cos(alpha_11) - U1*U2*abs(Y12)*cos(theta_1 - alpha_12);

fprintf('P1 = %0.4f + %0.4f*sin(theta_1 - %0.4f)\n', (U1^2)*abs(Y11)*sin(alpha_11), U1*U2*abs(Y12), alpha_12*toDeg)
fprintf('Q1 = %0.4f - %0.4f*cos(theta_1 - %0.4f)\n', (U1^2)*abs(Y11)*cos(alpha_11), U1*U2*abs(Y12), alpha_12*toDeg)

% Ъглови х-ки на системата
P2_theta = -(U2^2)*abs(Y22)*sin(alpha_22) - U2*U1*abs(Y21)*sin(-theta_1 - alpha_21);
Q2_theta = -(U2^2)*abs(Y22)*cos(alpha_22) + U2*U1*abs(Y21)*cos(-theta_1 - alpha_21);

fprintf('P2 = -%0.4f - %0.4f*sin(-theta_1 - %0.4f)\n', (U2^2)*abs(Y22)*sin(alpha_22), U2*U1*abs(Y21), alpha_21*toDeg)
fprintf('Q2 = -%0.4f + %0.4f*cos(-theta_1 - %0.4f)\n', (U2^2)*abs(Y22)*cos(alpha_22), U2*U1*abs(Y21), alpha_21*toDeg)

% Графики
figure()
subplot(1,2,1)
plot(theta_1*toDeg, P1_theta, 'k-', 'LineWidth', 1.5)
hold on; box on; grid on;
plot(theta_1*toDeg, P2_theta, 'k-.', 'LineWidth', 1.5)
xlabel '\theta_1 [deg]'
ylabel 'P [pu]'
legend({'P_1(\theta_1)', 'P_2(\theta_1)'})
xlim([0 180])

subplot(1,2,2)
plot(theta_1*toDeg, Q1_theta, 'k-', 'LineWidth', 1.5)
hold on; box on; grid on;
plot(theta_1*toDeg, Q2_theta, 'k-.', 'LineWidth', 1.5)
xlabel '\theta_1 [deg]'
ylabel 'Q [pu]'
legend({'Q_1(\theta_1)', 'Q_2(\theta_1)'})
xlim([0 180])

%% b) Установен режим на мрежата
fprintf('\nb) Установен режим на мрежата\n')

% Ъгъл на генераторното напрежение
theta_10 = asin( (P1 - (U1^2)*abs(Y11)*sin(alpha_11))/(U1*U2*abs(Y12)) ) + alpha_12;
fprintf('theta_1 = %0.4f deg\n', theta_10 * toDeg)

% Мощности на генератора и системата
P1 = (U1^2)*abs(Y11)*sin(alpha_11) + U1*U2*abs(Y12)*sin(theta_10 - alpha_12);
Q1 = (U1^2)*abs(Y11)*cos(alpha_11) - U1*U2*abs(Y12)*cos(theta_10 - alpha_12);
fprintf('P1 = %0.4f pu\n', P1)
fprintf('Q1 = %0.4f pu\n', Q1)

P2 = -(U2^2)*abs(Y22)*sin(alpha_22) - U2*U1*abs(Y21)*sin(-theta_10 - alpha_21);
Q2 = -(U2^2)*abs(Y22)*cos(alpha_22) + U2*U1*abs(Y21)*cos(-theta_10 - alpha_21);
fprintf('P2 = %0.4f pu\n', P2)
fprintf('Q2 = %0.4f pu\n', Q2)

% Консумация в мрежата (вкл. загуби)
fprintf('Консумация в мрежата включително загуби:\n')
fprintf('  * Aктивна мощност:   P1-P2 = %0.4f pu \n', P1-P2)
fprintf('  * Реактивна мощност: Q1-Q2 = %0.4f pu \n', Q1-Q2)

%% c) Установен режим на генератора
fprintf('\nc) Установен режим на генератора\n')

fprintf('\nОт установения режим на мрежата е известно:\n')
Ug0 = U1*exp(1i*theta_10);
fprintf('Ug0 = (%0.4f /_%0.4f deg) pu\n', abs(Ug0), angle(Ug0)*toDeg)

Pg0 = P1;
Qg0 = Q1;
fprintf('Pg0 = %0.4f pu;   Qg0 = %0.4f pu\n', P1, Q1)

fprintf('\nИзчисляват се:\n')
Ig0 = (Pg0 - 1i*Qg0) / conj(Ug0);
fprintf('Ig0 = I_R0 + jI_I0 = %0.4f %+0.4fi pu \n', real(Ig0), imag(Ig0))

delta0 = angle(Ug0 + (Ra+1i*Xq)*Ig0);

fprintf('delta0 = %0.4f deg\n', delta0 * toDeg)

fprintf('omega_r0 = 1 pu\n')

Id0 = real(Ig0)*sin(delta0) - imag(Ig0)*cos(delta0);
Iq0 = real(Ig0)*cos(delta0) + imag(Ig0)*sin(delta0);
fprintf('Id0 = %0.4f pu  Iq0 = %0.4f pu \n', Id0, Iq0)

Ud0 = real(Ug0)*sin(delta0) - imag(Ug0)*cos(delta0);
Uq0 = real(Ug0)*cos(delta0) + imag(Ug0)*sin(delta0);
fprintf('Ud0 = %0.4f pu  Uq0 = %0.4f pu \n', Ud0, Uq0)

Ed0_prim = Ud0 + Ra*Id0 - Xq_prim*Iq0;
Eq0_prim = Uq0 + Ra*Iq0 + Xd_prim*Id0;
fprintf('Ed0_prim = %0.4f pu  Eq0_prim = %0.4f pu \n', Ed0_prim, Eq0_prim)

Ed0_sec = Ud0 + Ra*Id0 - Xq_sec*Iq0;
Eq0_sec = Uq0 + Ra*Iq0 + Xd_sec*Id0;
fprintf('Ed0_sec  = %0.4f pu  Eq0_sec  = %0.4f pu \n', Ed0_sec, Eq0_sec)

psid0 =  Eq0_sec - Xd_sec*Id0;
psiq0 = -Ed0_sec - Xq_sec*Iq0;
fprintf('psid0    = %0.4f pu  psiq0    =%0.4f pu \n', psid0, psiq0)

Te0 = psid0*Iq0 - psiq0*Id0;
fprintf('Te0 = Pm0 = %0.4f pu \n', Te0)

EfD0 = K2*Eq0_prim - K1*Eq0_sec + K3*Id0;
fprintf('EfD0 = %0.4f pu \n', EfD0)

ER0_prim =  Ed0_prim*sin(delta0) + Eq0_prim*cos(delta0);
EI0_prim = -Ed0_prim*cos(delta0) + Eq0_prim*sin(delta0);
fprintf('ER0_prim = %0.4f pu  EI0_prim = %0.4f pu \n', ER0_prim, EI0_prim)

ER0_sec =  Ed0_sec*sin(delta0) + Eq0_sec*cos(delta0);
EI0_sec = -Ed0_sec*cos(delta0) + Eq0_sec*sin(delta0);
fprintf('ER0_sec  = %0.4f pu  EI0_sec  = %0.4f pu \n', ER0_sec, EI0_sec)

EQ0 = Eq0_prim + (Xq-Xd_prim)*Id0;
fprintf('EQ0  = %0.4f pu \n', EQ0)

% Графики
fsize = 10;
figure()
subplot(1,2,1)
plot(theta_1*toDeg, P1_theta, 'k-', 'LineWidth', 1.5)
hold on; box on; grid on;
plot([theta_10, theta_10]*toDeg,[0, P1],'k--o', 'LineWidth', 1.5, 'MarkerSize', 3 )
plot([0,        theta_10]*toDeg,[P1,P1],'k--o', 'LineWidth', 1.5, 'MarkerSize', 3 )
text(theta_10*toDeg+5, 0, '\theta_{0}','FontSize',fsize, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
text(-3, P1, 'P_{G0}','FontSize',fsize, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle')
xlabel '\theta_1 [deg]'
ylabel 'P_1 [pu]'
xlim([0 180])

subplot(1,2,2)
plot(theta_1*toDeg, Q1_theta, 'k-', 'LineWidth', 1.5)
hold on; box on; grid on;
plot([theta_10, theta_10]*toDeg,[0, Q1],'k--o', 'LineWidth', 1.5, 'MarkerSize', 3 )
plot([0,        theta_10]*toDeg,[Q1,Q1],'k--o', 'LineWidth', 1.5, 'MarkerSize', 3 )
text(theta_10*toDeg+5, 0, '\theta_{0}','FontSize',fsize, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
text(-3, Q1, 'Q_{G0}','FontSize',fsize, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle')
xlabel '\theta_1 [deg]'
ylabel 'Q_1 [pu]'
xlim([0 180])

%% d) Векторна диаграма на генератора
fprintf('\nc) Установен режим на генератора\n')

figure()
of = 0.1;
% Плотиране на R-I координатите
len = 2.5;
quiver(0,0, len, 0, 0, 'k-', 'LineWidth', 0.5)
hold on; grid on; box on;
quiver(0,0, 0, len, 0, 'k-', 'LineWidth', 0.5)
text(len+of, 0, 'R','FontSize',fsize, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle')
text(0, len+of, 'I','FontSize',fsize, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom')
axis equal

% Плотиране на d-q координатната система
xd = len*sin(delta0);
yd =-len*cos(delta0);
xq = len*cos(delta0);
yq = len*sin(delta0);
quiver(0,0, xd, yd, 0, 'k-', 'LineWidth', 0.5)
quiver(0,0, xq, yq, 0, 'k-', 'LineWidth', 0.5)
text(xd+of, yd, 'd','FontSize',fsize, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle')
text(xq, yq+of, 'q','FontSize',fsize, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom')
axis equal

% Плотиране на Ug и Ig
quiver(0,0, real(Ug0), imag(Ug0), 0, 'k-', 'LineWidth', 1.5)
text(real(Ug0)+of, imag(Ug0), 'U_{G0}','FontSize',fsize, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle')

quiver(0,0, real(Ig0), imag(Ig0), 0, 'k-', 'LineWidth', 1.5)
text(real(Ig0)+of, imag(Ig0), 'I_{G0}','FontSize',fsize, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom')

% Плотиране на E' и E"
quiver(0,0, ER0_prim, EI0_prim, 0, 'k-', 'LineWidth', 1.5)
text(ER0_prim+of, EI0_prim, 'E^{''}_0','FontSize',fsize, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom')

quiver(0,0, ER0_sec, EI0_sec, 0, 'k-', 'LineWidth', 1.5)
text(ER0_sec+of, EI0_sec, 'E^{''''}_0','FontSize',fsize, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle')

% Плотиране на EQ
quiver(0,0, EQ0*cos(delta0), EQ0*sin(delta0), 0, 'k-', 'LineWidth', 1.5)
text(EQ0*cos(delta0)+of, EQ0*sin(delta0), 'E_{Q0}','FontSize',fsize, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')

xlim([-0.5, 3])
ylim([-1.5, 3])

% Плотиране на ъгъл delta0
delta = linspace(0, delta0, 30);
R = 1.7; % радиус на дъгичката
arc_x = R*cos(delta);
arc_y = R*sin(delta);
plot(arc_x, arc_y, 'k-', 'LineWidth', 0.5)
text(arc_x(15)+of, arc_y(15), sprintf('%s = %0.3f%s', '\delta_0', delta0*toDeg, char(176)),...
    'FontSize',fsize, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle')