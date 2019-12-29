clear all; clc; close all;
format shortg;

datafile = 'raw_gps_log2';
data = load(datafile);

idx = data(:, 1) != 0.0;
data = data(idx, :);

printf('lattitude,longitude,speed,heading\n');

for i = 1:length(data)
    lat = data(i, 1);
    lon = data(i, 2);
    spd_kph = data(i, 3);
    hdg = data(i, 4);

    % degrees = sign(lat)*floor(abs(lat) / 100);
    % minutes = (lat - 100*degrees);
    % decimal = minutes / 60;
    % lat = degrees + decimal;

    % degrees = sign(lon)*floor(abs(lon) / 100);
    % minutes = (lon - 100*degrees);
    % decimal = minutes / 60;
    % lon = degrees + decimal;

    lat = decimalminutes2decimal(lat);
    lon = decimalminutes2decimal(lon);

    printf('%f,%f,%f,%f\n', lat, lon, spd_kph, hdg);
end
