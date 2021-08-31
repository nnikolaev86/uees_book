%% ������ 2: ����������� �� ������ �������������� � ��������� ����� �� ����������� ���
clear all; clc;
fprintf('������ 2: ����������� �� ������ �������������� � ��������� ����� �� ����������� ���\n')

%% ������ �����

%��������� ������������:
Y11  = 0.1882 -2.6990i ;
Y22  = 0.1802 -2.8672i ;
%������� ������������:
Y12  = 0.1407 -2.5663i ;
Y21  = Y12;

% �� ������������ �����
U1 = 1.05;
P1 = 0.7;
% �� ���������
U2 = 1.0;
theta2 = 0; % ������ �� ��������� � ���� ������ �� �� �������� ������ � �����������

% ����. �� ������������� �� RAD � DEG
toDeg = 180/pi; 

% ��������� �� ����������
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

%% �) ������ ��������������
fprintf('\na) ������ ��������������\n')

% ��������� ���� �� �������������� (������ �� � RAD)
alpha_11 = pi/2 + angle(Y11);
alpha_22 = pi/2 + angle(Y22);
alpha_12 = pi/2 + angle(Y12);
alpha_21 = pi/2 + angle(Y21);

fprintf('alpha_11 = %0.4f deg\n', alpha_11 * toDeg)
fprintf('alpha_22 = %0.4f deg\n', alpha_22 * toDeg)
fprintf('alpha_12 = %0.4f deg\n', alpha_12 * toDeg)
fprintf('alpha_21 = %0.4f deg\n\n', alpha_21 * toDeg)

% �������� �� ��������� �� ���� theta_1
theta_1 = linspace(0, pi, 100); % ���������� �� 100 ����� ����� 0 � pi/2 RAD (0 �� 90 DEG)

% ������ �-�� �� ����������
P1_theta = (U1^2)*abs(Y11)*sin(alpha_11) + U1*U2*abs(Y12)*sin(theta_1 - alpha_12);
Q1_theta = (U1^2)*abs(Y11)*cos(alpha_11) - U1*U2*abs(Y12)*cos(theta_1 - alpha_12);

fprintf('P1 = %0.4f + %0.4f*sin(theta_1 - %0.4f)\n', (U1^2)*abs(Y11)*sin(alpha_11), U1*U2*abs(Y12), alpha_12*toDeg)
fprintf('Q1 = %0.4f - %0.4f*cos(theta_1 - %0.4f)\n', (U1^2)*abs(Y11)*cos(alpha_11), U1*U2*abs(Y12), alpha_12*toDeg)

% ������ �-�� �� ���������
P2_theta = -(U2^2)*abs(Y22)*sin(alpha_22) - U2*U1*abs(Y21)*sin(-theta_1 - alpha_21);
Q2_theta = -(U2^2)*abs(Y22)*cos(alpha_22) + U2*U1*abs(Y21)*cos(-theta_1 - alpha_21);

fprintf('P2 = -%0.4f - %0.4f*sin(-theta_1 - %0.4f)\n', (U2^2)*abs(Y22)*sin(alpha_22), U2*U1*abs(Y21), alpha_21*toDeg)
fprintf('Q2 = -%0.4f + %0.4f*cos(-theta_1 - %0.4f)\n', (U2^2)*abs(Y22)*cos(alpha_22), U2*U1*abs(Y21), alpha_21*toDeg)

% �������
figure()
plot(theta_1*toDeg, P1_theta, 'k-', 'LineWidth', 1.5)
hold on; box on; grid on;
plot(theta_1*toDeg, P2_theta, 'k-.', 'LineWidth', 1.5)
xlabel '\theta_1 [deg]'
ylabel 'P [pu]'
legend({'P_1(\theta_1)', 'P_2(\theta_1)'})

figure()
plot(theta_1*toDeg, Q1_theta, 'k-', 'LineWidth', 1.5)
hold on; box on; grid on;
plot(theta_1*toDeg, Q2_theta, 'k-.', 'LineWidth', 1.5)
xlabel '\theta_1 [deg]'
ylabel 'Q [pu]'
legend({'Q_1(\theta_1)', 'Q_2(\theta_1)'})

%% b) ��������� ����� �� �������
fprintf('\nb) ��������� ����� �� �������\n')

% ���� �� ������������� ����������
theta_1 = asin( (P1 - (U1^2)*abs(Y11)*sin(alpha_11))/(U1*U2*abs(Y12)) ) + alpha_12;
fprintf('theta_1 = %0.4f deg\n', theta_1 * toDeg)

% �������� �� ���������� � ���������
P1 = (U1^2)*abs(Y11)*sin(alpha_11) + U1*U2*abs(Y12)*sin(theta_1 - alpha_12);
Q1 = (U1^2)*abs(Y11)*cos(alpha_11) - U1*U2*abs(Y12)*cos(theta_1 - alpha_12);
fprintf('P1 = %0.4f pu\n', P1)
fprintf('Q1 = %0.4f pu\n', Q1)

P2 = -(U2^2)*abs(Y22)*sin(alpha_22) - U2*U1*abs(Y21)*sin(-theta_1 - alpha_21);
Q2 = -(U2^2)*abs(Y22)*cos(alpha_22) + U2*U1*abs(Y21)*cos(-theta_1 - alpha_21);
fprintf('P2 = %0.4f pu\n', P2)
fprintf('Q2 = %0.4f pu\n', Q2)

% ���������� � ������� (���. ������)
fprintf('���������� � ������� ����������� ������:\n')
fprintf('  * A������ �������:   P1-P2 = %0.4f pu \n', P1-P2)
fprintf('  * ��������� �������: Q1-Q2 = %0.4f pu \n', Q1-Q2)

%% c) ��������� ����� �� ����������
fprintf('\nc) ��������� ����� �� ����������\n')

fprintf('\n�� ����������� ����� �� ������� � ��������:\n')
Ug0 = U1*exp(1i*theta_1);
fprintf('Ug0 = (%0.4f /_%0.4f deg) pu\n', abs(Ug0), angle(Ug0)*toDeg)

Pg0 = P1;
Qg0 = Q1;
fprintf('Pg0 = %0.4f pu;   Qg0 = %0.4f pu\n', P1, Q1)

fprintf('\n���������� ��:\n')
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

ER0_sec =  Ed0_sec*sin(delta0) + Eq0_sec*cos(delta0);
EI0_sec = -Ed0_sec*cos(delta0) + Eq0_sec*sin(delta0);
fprintf('ER0_sec  = %0.4f pu  EI0_sec  = %0.4f pu \n', ER0_sec, EI0_sec)















