clear all, close all, clc

% Transmitting coil
V = 100; % filament current [A]
turns = 410; % Numer of turn in the coil []
radius = 1; % Radius of transmitting coil [mm?]
wireGauge = 6;

% Receiving coil
turns_car = 320; % Numer of turn in the coil []
radius_car = 1.5; % Radius of transmitting coil [mm?]
wireGauge_car = 8;

% Coil orientation
height = .15; % Height of the receiving coil above the transmitting coil MUST BE ON THE GRID
spacing = .5; % distance betewen centers of coils
velocity = 30; % velocity of car [m/s]

%% MAIN

totalScenarios = 2;
data = zeros(totalScenarios, 11); % Initialize data matrix
format shortG % print data with the desired detail

%totalCharge = DWPT(V, wireGauge, turns, radius, wireGauge_car, turns_car, radius_car, height, spacing, velocity, 1);
%data(1,:) = [V wireGauge turns radius wireGauge_car turns_car radius_car height spacing velocity totalCharge];

totalCharge = DWPT(600, wireGauge, 510, radius, wireGauge_car, turns_car, radius_car, height, spacing, velocity, 2);
data(2,:) = [600 wireGauge 510 radius wireGauge_car turns_car radius_car height spacing velocity totalCharge];

%writematrix(data,'dataComp/data.out')