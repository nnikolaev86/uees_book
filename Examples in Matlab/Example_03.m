%% ПРИМЕР 3: Устойчивост при малки смущения на модел от 2-ри ред при вариращо Eq_prim
clear all; clc; close all;
fprintf('ПРИМЕР 3: Устойчивост при малки смущения на модел от 2-ри ред при вариращо Eq_prim\n')

%% Входни данни

% Режимни параметри на генератора и системата
U1     = 1.02;
P1     = 0.2;

U2     = 1.0;
theta2 = 0; % ъгълът на системата е нула затова ще се пропуска надолу в уравненията

% Номинална кръгова честота
w0 = 2*pi*50; % rad/s

% Параметри на генератора за модел от 2-ри ред
Sn = 500; % MVA
H = 3.0;
Xd_prim  = 0.3;
Ra       = 0.0;

% Собствени и взаимни проводимости без включен генераторен импеданс
% (Това са изчислените в Example_01.m подточка с) проводимости)
Y11  = 0.1978 -2.5993i ;
Y22  = 0.1688 -2.4178i ;
Y12  = 0.1362 -2.7317i ;
Y21  = Y12;

alpha_11 = pi/2 + angle(Y11);
alpha_22 = pi/2 + angle(Y22);
alpha_12 = pi/2 + angle(Y12);
alpha_21 = pi/2 + angle(Y21);

% Собствени и взаимни проводимости с включен генераторен импеданс
% (изчислени са със скрипт Example_03_compute_admitances.m)
Y11g  = 0.0624 -1.4625i ;
Y22g  = 0.0853 -1.1603i ;
Y12g  = 0.0253 -1.5357i ;
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

%% b) Установен режим на генераторния модел
fprintf('\nb) Установен режим на генераторния модел\n')

delta0 = angle(Ug0 + (Ra + 1i*Xd_prim)*Ig0);
fprintf('delta0 = %0.4f deg\n', delta0 * toDeg)

Id0 = real(Ig0)*sin(delta0) - imag(Ig0)*cos(delta0);
Iq0 = real(Ig0)*cos(delta0) + imag(Ig0)*sin(delta0);
fprintf('Id0 = %0.4f pu  Iq0 = %0.4f pu \n', Id0, Iq0)

Ud0 = real(Ug0)*sin(delta0) - imag(Ug0)*cos(delta0);
Uq0 = real(Ug0)*cos(delta0) + imag(Ug0)*sin(delta0);
fprintf('Ud0 = %0.4f pu  Uq0 = %0.4f pu \n', Ud0, Uq0)

Eq0_prim = Uq0 + Ra*Iq0 + Xd_prim*Id0;
fprintf('Eq0_prim = %0.4f pu\n', Eq0_prim)

psid0 = Eq0_prim - Xd_prim*Id0;
psiq0 = -Xd_prim*Iq0;

Pm0 = psid0*Iq0 - psiq0*Id0;
fprintf('Pm0 = %0.4f pu\n', Pm0)

%% c) Линеаризиране на модела
fprintf('\nc) Линеаризиране на модела\n')

Ks = Eq0_prim*U2*abs(Y12g)*cos(delta0 - alpha_12g);
Kd = Pm0;

fprintf('Ks = %0.4f pu-в.момент/rad\n', Ks)
fprintf('Kd = %0.4f pu-в.момент/pu-скорост\n', Kd)

A = [        0,        w0;
     -Ks/(2*H), -Kd/(2*H)  ];

B = [     0;
     1/(2*H) ];

C = [1, 0;
     0, 1 ];
 
D = [0;
     0 ];

% Създаване на модел в пространство на състоянията
lin_model = ss(A, B, C, D);

% Задаване имената на входовете, изходите и променливите на състоянието
lin_model.StateName  = {'delta','w'};
lin_model.OutputName = {'delta','w'};
lin_model.InputName  = {'Pm'};

% Принтиране на модела
lin_model

%% d) Анализ на устойчивостта
fprintf('\nd) Анализ на устойчивостта\n')

% Корени
lambda = eig(A)

% Естествена недемпфирана честота
wn = sqrt(Ks*w0/(2*H));
fprintf('Естествена недемпфирана честота wn = %0.4f rad/s\n', wn)

% Декремент на затихване
xi = 0.5*Kd/sqrt(Ks*2*H*w0);
fprintf('Декремент на затихване xi = %0.4f \n', xi)

% Честота на колебанията
w = imag( lambda(1) );
fprintf('Честота на колебанията w = %0.4f rad/s = %0.4f Hz\n', w, w/(2*pi) )

%% Локус на корените
% Тук нареждаме корените за всички работни режими от A до Е
% тези с положителна имагинерна част
lambdas_pos = [ -0.0833 + 8.4256i;
                -0.0667 + 8.6243i;
                -0.0500 + 8.7727i;
                -0.0333 + 8.8751i;
                -0.0167 + 8.9342i;
                                ];
% и онези с отрицателната имагинерна част
lambdas_neg = [ -0.0833 - 8.4256i
                -0.0667 - 8.6243i;
                -0.0500 - 8.7727i;
                -0.0333 - 8.8751i;
                -0.0167 - 8.9342i;
                                ];
                            
cases = {'A','B','C','D','E'};
                            
% Графика на корените                            
plot(real(lambdas_pos), imag(lambdas_pos), 'k--', 'LineWidth', 0.2, 'Marker', 'x', 'MarkerSize', 7)
grid on; hold on
plot(real(lambdas_neg), imag(lambdas_neg), 'k--', 'LineWidth', 0.2, 'Marker', 'x', 'MarkerSize', 7)
xmin = -0.1;
xmax = 0.01;
ymin = -10;
ymax =  10;
ylim([ymin, ymax])
xlim([xmin, xmax])
xlabel 'Real   \alpha [Np/s]'
ylabel 'Imag   \omega [rad/s]'
plot([xmin, xmax], [0 0 ], 'k-', 'LineWidth', 0.1)
plot([0 0 ], [ymin, ymax], 'k-', 'LineWidth', 0.1)

for c=1:length(cases)
   text( real(lambdas_pos(c)) , imag(lambdas_pos(c)), cases(c),...
         'VerticalAlignment', 'top', 'HorizontalAlignment', 'center') 
     
     text( real(lambdas_neg(c)) , imag(lambdas_neg(c)), cases(c),...
         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center')     
end
