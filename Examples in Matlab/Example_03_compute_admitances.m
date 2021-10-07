% ������ �� ��������� �� ����������� � ������� ������������ � �������
% �������� �� ���������� ��� ��� �� ������ �� Example_03.m

clear all; clc
%% ������ �����
% ���������
Xd_prim  = 0.3;
Ra       = 0.0;
ZG       = Ra + 1i*Xd_prim;

% �������������
ZT = 0.15i; % pu

% �����
ZB = 11.11; % pu

% ��������������
ZW1 = 0.0225 + 0.225i;
YW1 = 0.11i;

ZW2 = 0.0225 + 0.225i;
YW2 = 0.11i;

ZW3 = 0.0225 + 0.225i;
YW3 = 0.11i;

ZW4 = 0.0225 + 0.225i;
YW4 = 0.11i;

%% ����� � ���������� ��������� �� ��������� ������������, ��� �������� �� ���������� ������������ �� ����������������
fprintf('\n����� � ���������� ��������� �� ��������� ������������, ��� ��������\n   ���������� ������������ �� ����������������\n')

% ��������� �� ��������� �� �������������� �� ��������� ��������
MW1 = [ (1/ZW1),        -(1/ZW1 + YW1/2);
    -(1/ZW1 + YW1/2), (1/ZW1)
    ];

MW2 = [ (1/ZW2),        -(1/ZW2 + YW2/2);
    -(1/ZW2 + YW2/2), (1/ZW2)
    ];

MW3 = [ (1/ZW3),        -(1/ZW3 + YW3/2);
    -(1/ZW3 + YW3/2), (1/ZW3)
    ];

MW4 = [ (1/ZW4),        -(1/ZW4 + YW4/2);
    -(1/ZW4 + YW4/2), (1/ZW4)
    ];

MT = [ 1/ZT, -1/ZT;
       -1/ZT,  1/ZT;
    ];

MG = [ 1/ZG, -1/ZG;
       -1/ZG,  1/ZG;
    ];

% ��������� �� ��������� �� �������������� �� �������
% � ������� 5�5, ������ ���� ������� � ���, ���� ��������� �����������
% ����� �� ������ �� ����� 5
Y = zeros(5,5);

% ��������� � ��� T
bus1 = 1;
bus2 = 3;
Y(bus1,bus1) = Y(bus1,bus1) + MT(1,1);
Y(bus2,bus2) = Y(bus2,bus2) + MT(2,2);
Y(bus1,bus2) = Y(bus1,bus2) + MT(1,2);
Y(bus2,bus1) = Y(bus2,bus1) + MT(2,1);

% ��������� � ��� W1
bus1 = 3;
bus2 = 4;
Y(bus1,bus1) = Y(bus1,bus1) + MW1(1,1);
Y(bus2,bus2) = Y(bus2,bus2) + MW1(2,2);
Y(bus1,bus2) = Y(bus1,bus2) + MW1(1,2);
Y(bus2,bus1) = Y(bus2,bus1) + MW1(2,1);

% ��������� � ��� W2
bus1 = 3;
bus2 = 4;
Y(bus1,bus1) = Y(bus1,bus1) + MW2(1,1);
Y(bus2,bus2) = Y(bus2,bus2) + MW2(2,2);
Y(bus1,bus2) = Y(bus1,bus2) + MW2(1,2);
Y(bus2,bus1) = Y(bus2,bus1) + MW2(2,1);

% ��������� � ��� W3
bus1 = 4;
bus2 = 2;
Y(bus1,bus1) = Y(bus1,bus1) + MW3(1,1);
Y(bus2,bus2) = Y(bus2,bus2) + MW3(2,2);
Y(bus1,bus2) = Y(bus1,bus2) + MW3(1,2);
Y(bus2,bus1) = Y(bus2,bus1) + MW3(2,1);

% ��������� � ��� W4
bus1 = 4;
bus2 = 2;
Y(bus1,bus1) = Y(bus1,bus1) + MW4(1,1);
Y(bus2,bus2) = Y(bus2,bus2) + MW4(2,2);
Y(bus1,bus2) = Y(bus1,bus2) + MW4(1,2);
Y(bus2,bus1) = Y(bus2,bus1) + MW4(2,1);

% ��������� � ��� B
bus1 = 3;
Y(bus1,bus1) = Y(bus1,bus1) + 1/ZB;

% ��������� � ��� �
bus1 = 1;
bus2 = 5;
Y(bus1,bus1) = Y(bus1,bus1) + MG(1,1);
Y(bus2,bus2) = Y(bus2,bus2) + MG(2,2);
Y(bus1,bus2) = Y(bus1,bus2) + MG(1,2);
Y(bus2,bus1) = Y(bus2,bus1) + MG(2,1);

fprintf('������� � �������������� �� �������')
Y

% ��������� ������� ���������� � � ��� �����������
% ������������ ������� � �������������� �� �������,
% � ����� ���������� ���� ����� 1 � 2
bus_G = [2 5];
bus_B = [1 3 4];

YGG = Y(bus_G, bus_G);
YGB = Y(bus_G, bus_B);
YBG = Y(bus_B, bus_G);
YBB = Y(bus_B, bus_B);

fprintf('������������ ������� � ��������� ������������ �')
YG = YGG - YGB*inv(YBB)*YBG

% �� ������������ ������� ��������� ����������� � �������
% ������������ �� ����� ����������� �����
% ��������� ������������
fprintf('��������� ������������:\n')
Y11g = YG(2,2);
fprintf('Y11g  = %0.4f %+0.4fi \n', real(Y11g), imag(Y11g))

Y22g = YG(1,1);
fprintf('Y22g  = %0.4f %+0.4fi \n', real(Y22g), imag(Y22g))

% ������� ������������
fprintf('������� ������������:\n')
Y12g = -YG(1,2);

fprintf('Y12g  = %0.4f %+0.4fi \n', real(Y12g), imag(Y12g))
Y21g = -YG(2,1);
fprintf('Y21g  = %0.4f %+0.4fi \n', real(Y21g), imag(Y21g))