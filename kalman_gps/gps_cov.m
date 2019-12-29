clear all; clc; close all;

deg2rad = @(d) pi*d/180;
rad2deg = @(r) 180*r/pi;

data = load('TrackWaypoint');
Re = 6.3781e6; % Radius of earth, m

lat  = deg2rad(data(:, 1));
lon  = deg2rad(data(:, 2));
alt  = data(:, 3);
t    = data(:, 4) - data(1, 4);
vel  = data(:, 5) * (1000/3600); % km/h to m/s
hdg  = data(:, 6);

state_data = [lat, lon, vel];
DIMS = size(state_data, 2);

% P = cov(state_data);
P = eye(DIMS);
Q = 0.1 * eye(DIMS);

H = eye(DIMS);

pred_output = zeros(length(lat), DIMS);

xhat = state_data(1, :)';
pred_output(1, :) = xhat;

for i = 2:length(lat)
    dt = t(i) - t(i-1);

    ctheta = cosd(hdg(i));
    stheta = sind(hdg(i));

    % Predict, with velocity and heading
    F = [1    0   ctheta*dt/Re;  % [ lat ]
         0    1   stheta*dt/Re;  % [ lon ]
         0    0       1       ]; % [ vel ]

    xhat_next = F * xhat;
    P = F * P * F' + Q;

    % Expected measurement
    mu_expect    = H * xhat_next;
    Sigma_expect = H * P * H';

    % actual measurement
    mu_actual = state_data(i, :)';
    Sigma_actual = cov(state_data);

    % kalman gain
    K = Sigma_expect * (Sigma_expect + Sigma_actual)^-1;

    % update
    xhat_next = xhat_next + K * (mu_actual - mu_expect);
    P = P - K * H * P;
    pred_output(i, :) = xhat_next;

    xhat = xhat_next;
end

pred_lat = pred_output(:, 1);
pred_lon = pred_output(:, 2);

figure; hold on; grid on;
title('Kalman vs sensor, with velocity data');
xlabel('Longitude');
xlabel('Lattitude');
plot(rad2deg(lon),      rad2deg(lat),      'b', 'LineWidth', 2)
plot(rad2deg(pred_lon), rad2deg(pred_lat), 'r', 'LineWidth', 2)


pred_output = zeros(length(lat), DIMS);

xhat = state_data(1, :)';
pred_output(1, :) = xhat;

for i = 2:length(lat)
    dt = t(i) - t(i-1);

    % Predict without considering velocity or heading
    F = eye(DIMS);

    xhat_next = F * xhat;
    P = F * P * F' + Q;

    % Expected measurement
    mu_expect    = H * xhat_next;
    Sigma_expect = H * P * H';

    % actual measurement
    mu_actual = state_data(i, :)';
    Sigma_actual = P;

    % kalman gain
    K = Sigma_expect * (Sigma_expect + Sigma_actual)^-1;

    % update
    xhat_next = xhat_next + K * (mu_actual - mu_expect);
    P = P - K * H * P;
    pred_output(i, :) = xhat_next;

    xhat = xhat_next;
end

pred_lat = pred_output(:, 1);
pred_lon = pred_output(:, 2);

figure; hold on; grid on;
width = 1; % line width
title('Kalman vs sensor, no velocity data');
xlabel('Longitude');
xlabel('Lattitude');
plot(rad2deg(lon),      rad2deg(lat),      'b', 'LineWidth', 1)
plot(rad2deg(pred_lon), rad2deg(pred_lat), 'r', 'LineWidth', 1)
try
    c_data = load('c_kalman_output');
    c_lat = c_data(:, 1);
    c_lon = c_data(:, 2);
    plot(c_lon, c_lat, 'g', 'LineWidth', 1)
catch
    disp('Could not find C program output');
end
print('no_vel_kalman', '-dpng');

