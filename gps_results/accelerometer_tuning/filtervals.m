clear all; clc; close all;
format shortg;

stationary = load('stationarydump.txt');

t  = stationary(:, 1);
ax = stationary(:, 2);
ay = stationary(:, 3);
az = stationary(:, 4);
gx = stationary(:, 5);
gy = stationary(:, 6);
gz = stationary(:, 7);

fprintf('AX: add %3.8f to bring to zero, sigma %3.8f\n', -mean(ax), std(ax));
fprintf('AY: add %3.8f to bring to zero, sigma %3.8f\n', -mean(ay), std(ay));
fprintf('AZ: add %3.8f to bring to 0.5,  sigma %3.8f\n', 0.5-mean(az), std(az));

fprintf('GX: add %3.8f to bring to zero, sigma %3.8f\n', -mean(gx), std(gx));
fprintf('GY: add %3.8f to bring to zero, sigma %3.8f\n', -mean(gy), std(gy));
fprintf('GZ: add %3.8f to bring to zero, sigma %3.8f\n', -mean(gz), std(gz));


moving = load('movingdump.txt');
t  = moving(:, 1);
ax = moving(:, 2);
ay = moving(:, 3);
az = moving(:, 4);
gx = moving(:, 5);
gy = moving(:, 6);
gz = moving(:, 7);


figure; hold on; grid on;
plot(t, ax, 'linewidth', 2, 'r');
plot(t, ay, 'linewidth', 2, 'g');
plot(t, az, 'linewidth', 2, 'b');
% pause

figure; hold on; grid on;
plot(t, gx, 'linewidth', 2, 'r');
plot(t, gy, 'linewidth', 2, 'g');
plot(t, gz, 'linewidth', 2, 'b');

vector = [ax, ay, az, gx, gy, gz];
disp('moving covariance');
disp(cov(vector));

