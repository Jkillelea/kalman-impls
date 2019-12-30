clear all; clc; close all;

datafile = './gps_results/accelerometer_tuning/accelerometerdump';

data = load(datafile);
datasize = size(data, 1);

t = data(:, 1);

% AX: add -0.03918230 to bring to zero, sigma 0.00171112
% AY: add  0.00270217 to bring to zero, sigma 0.00151103
% AZ: add  0.03682770 to bring to 0.5,  sigma 0.00258348
% GX: add -0.00972417 to bring to zero, sigma 0.00050139
% GY: add  0.00645349 to bring to zero, sigma 0.00041639
% GZ: add  0.00423421 to bring to zero, sigma 0.00041423
a = 9.81 * (data(:, 2:4) + repmat([-0.03918230, 0.00270217, 0.03682770], [datasize, 1]));
g = data(:, 5:7) + 250*repmat([-0.00972417, 0.00645349, 0.00423421], [datasize, 1]);


DIMS = 12;
x       = zeros(DIMS, 1);
mu_e    = zeros(DIMS, 1);
mu_a    = zeros(DIMS, 1);
P       = 0.1*eye(DIMS, DIMS);
Q       = 0.1*eye(DIMS, DIMS);
H       = eye(DIMS, DIMS);
Sigma_e = zeros(DIMS, DIMS);
Sigma_a = 0.001*eye(DIMS, DIMS);
K       = zeros(DIMS, DIMS);
F       = eye(DIMS, DIMS);

dt = 0.1;
F(1, 7)  = dt;
F(2, 8)  = dt;
F(3, 9)  = dt;
F(4, 10) = dt;
F(5, 11) = dt;
F(6, 12) = dt;

for i = 1:6
    H(i, i) = 0;
end

xdata = zeros(datasize, DIMS);

for i = 1:datasize
    % fprintf('%d/%d\n', i, datasize);

    phi   = x(4);
    theta = x(5);
    psi   = x(6);

    DCM1 = [
        1   0           0;
        0   cosd(phi)   sind(phi);
        0  -sind(phi)   cosd(phi);
    ];

    DCM2 = [
        cosd(theta)  0  -sind(theta);
        0            1            0;
        sind(theta)  0   cosd(theta);
    ];

    DCM3 = [
         cosd(psi)  sind(psi)  0;
        -sind(psi)  cosd(psi)  0;
         0          0          1;
    ];

    DCM = DCM3 * DCM2 * DCM1;

    x = F * x;
    P = F * (P * F') + Q;

    mu_e    = H*x;
    Sigma_e = H * (P * H');

    mu_a(7:9)   = DCM * a(i, :)';
    mu_a(10:12) = DCM * g(i, :)';

    K = P * H' * ((Sigma_e - Sigma_a)^-1);

    x = x + K * (mu_a - mu_e);
    P = P - K * (H * P);

    xdata(i, :) = x;
end

figure; hold on; grid on;
title('ax');
plot(t, a(:, 1),      'r');
plot(t, xdata(:, 7), 'b');
plot(t, xdata(:, 1), 'b');

figure; hold on; grid on;
title('ay');
plot(t, a(:, 2),      'r');
plot(t, xdata(:, 8), 'b');
plot(t, xdata(:, 2), 'b');

figure; hold on; grid on;
title('az');
plot(t,     a(:, 3), 'r');
plot(t, xdata(:, 9), 'b');
plot(t, xdata(:, 3), 'b');

figure; hold on; grid on;
title('gx');
plot(t, g(:, 1),      'r');
plot(t, xdata(:, 10), 'b');
plot(t, xdata(:, 4), 'b');

figure; hold on; grid on;
title('gy');
plot(t, g(:, 2),      'r');
plot(t, xdata(:, 11), 'b');
plot(t, xdata(:, 5), 'b');

figure; hold on; grid on;
title('gz');
plot(t, g(:, 3),      'r');
plot(t, xdata(:, 12), 'b');
plot(t, xdata(:, 6), 'b');




