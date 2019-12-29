clear all; clc; close all;

deg2rad = @(d) pi*d/180;
rad2deg = @(r) 180*r/pi;

data = load('cleaned_raw_gps_octavefmt');
Re = 6.3781e6; % Radius of earth, meters

lat  = deg2rad(data(:, 1));
lon  = deg2rad(data(:, 2));
vel  = data(:, 3) * (1000/3600); % km/h to m/s
hdg  = data(:, 4);

% mu_actual at each time interval
state_data = [lat, lon, vel];
DIMS = size(state_data, 2);

% Predicted measurement covariance
% P = cov(state_data)
P = 99999.9*eye(DIMS);

% Increase in variance per iteration
Q = 0.1 * eye(DIMS);

% expect 1:1 between what's measured and what's reported
H = eye(DIMS);

pred_output = zeros(length(lat), DIMS);

xhat = state_data(1, :)';
pred_output(1, :) = xhat;

fid = fopen('gps_kalman_result', 'w');
fprintf(fid, 'lattitude,longitude\n');
fprintf(fid, '%f,%f\n', rad2deg(xhat(1)), rad2deg(xhat(2)));

for i = 2:length(lat)
    dt = 1;

    % Predict, with velocity and heading update latitude and longitude
    F = [1    0   cosd(hdg(i))*dt/Re; % [ lat ]
         0    1   sind(hdg(i))*dt/Re; % [ lon ]
         0    0   1                ]; % [ vel ]

    xhat_next = F * xhat;
    P = F * P * F' + Q;

    % Expected measurement
    mu_expect    = H * xhat_next;
    Sigma_expect = H * P * H';

    % actual measurement
    mu_actual = state_data(i, :)';
    % Observed measurement covariance
    % TODO: make this a function of PDOP
    Sigma_actual = 0.1*eye(3);
    Sigma_actual(3, 3) = P(3, 3);

    % kalman gain
    K = Sigma_expect * (Sigma_expect + Sigma_actual)^-1;

    % update
    xhat_next = xhat_next + K * (mu_actual - mu_expect);
    P = P - K * H * P;
    pred_output(i, :) = xhat_next;

    xhat = xhat_next;

    fprintf(fid, '%f,%f\n', rad2deg(xhat(1)), rad2deg(xhat(2)));
end

fclose(fid);
