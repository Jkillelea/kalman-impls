% experiment with noise tolerance and sample rates for a kalman filter
clear all;
close all;
clc;

% octave required
pkg load odepkg

% truth data is very high frequency w/ no noise
dt = 1;
tmax = 100;
t = 0:dt:tmax;

% velocity function
v = @(t) 2*sin(0.1*t) + 0.1;

[t, x_pos] = ode45(v, t, 0);

vel = v(t);

state = [x_pos, vel];

% samples at lower frequency and add noise
pos_err  = 10;
pos_bias = 2;
vel_err  = 1;
vel_bias = 0.1;

pos_unc = pos_err + pos_bias;
vel_unc = vel_err + vel_bias;

x_sample = x_pos ...
         + pos_err*(rand(length(t), 1) - 0.5) ...
         + pos_bias*(rand() - 0.5);

vel_sample = vel ...
           + vel_err*(rand(length(t), 1) - 0.5) ...
           + vel_bias*(rand() - 0.5);

sample_state = [x_sample, vel_sample];

% Calculate kalman filter for each sample point
% forward propagation matrix
F = [1 dt;
     0 1];
% Covariance matrix (zero assuming initial state is known exactly)
P = zeros(2);
% External uncertainty sources
Q = 0.1*ones(2);
% (initial) state vector
pred_state = [0; 0];
% observation matrix
H = [1, 0;
     0, 1];
% Measurement Uncertainty Covariances
R = [pos_unc, 0;
     0,       vel_unc];

predicted_states = zeros(length(t), 2);
for i = 2:length(t)
    dt = t(i) - t(i-1);

    % Predicted next state
    pred_state = F * pred_state;
    % Next covariance (adds uncertainty)
    P = F * P * F' + Q;

    % Predicted measurement
    mu0 = H * pred_state;
    % Predicted measurement covariance
    Sigma0 = H * P * H';

    % Observed measurement
    z   = sample_state(i, :)';
    mu1 = z;
    % Observed measurement uncertainty covariance
    Sigma1 = R;

    % Kalman gain
    K = Sigma0 * (Sigma0 + Sigma1)^-1;

    % update state prediction based on expected vs actual measurement
    pred_state = pred_state + K * (mu1 - mu0);
    P = P - K * H * P;

    predicted_states(i, :) = pred_state;
end

% plot x vs time
figure; hold on; grid on;
title('position');
plot(t, x_pos, 'linewidth', 2, 'displayname', 'true x')
plot(t, x_sample, 'displayname', 'samples')
plot(t, predicted_states(:, 1), 'linewidth', 2, 'displayname', 'kalman x')
legend('show')

% plot v vs time
figure; hold on; grid on;
title('velocity');
plot(   t, vel, 'linewidth', 2, 'displayname', 'true v')
plot(t, vel_sample, 'displayname', 'samples')
plot(   t, predicted_states(:, 2), 'linewidth', 2, 'displayname', 'kalman x')
legend('show')

% plot x vs v
figure; hold on; grid on;
title('pos vs vel');
plot(vel, x_pos, 'linewidth', 2, 'displayname', 'x\_pos vs vel')
plot(vel_sample, x_sample, 'displayname', 'x\_sample vs vel\_sample')
plot(predicted_states(:, 2), predicted_states(:, 1), 'linewidth', 2, ...
                                        'displayname', 'kalman v vs x')
legend('show')

