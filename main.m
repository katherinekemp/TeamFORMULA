clear all, close all, clc

%% FORMULA variables

% Transmitting coil
V = 10; % Voltage of solar panels [V]
wireGauge = 8; % Diameter of wire [m]
turns = 10; % Numer of turn in the coil []
radius = 10; % Radius of transmitting coil [m]

% Receiving coil
wireGauge_car = 12; % Diameter of wire [m]
turns_car = 10; % Numer of turn in the coil []
radius_car = 10; % Radius of transmitting coil [m]

% Coil orientation
height = 2; % Height of the receiving coil above the transmitting coil [m]
spacing = 3 * radius; % distance betewen centers of [m]
velocity = 1; % velocity of car [m/s]

%% Cost Analysis Variables


%% MAIN

totalScenarios = 2;
data = zeros(totalScenarios, 11); % Initialize data matrix
format shortG % print data with the desired detail

%{
for i = 1:totalScenarios
    scenarioID = 2;
    totalCharge = DWPTeff(V, wireGauge, turns, radius, wireGauge_car, turns_car, radius_car, height, spacing, velocity, scenarioID);
    data(i,:) = [V wireGauge turns radius wireGauge_car turns_car radius_car height spacing velocity totalCharge];
end
%}
totalCharge = DWPT(V, wireGauge, turns, radius, wireGauge_car, turns_car, radius_car, height, spacing, velocity, 1);
data(1,:) = [V wireGauge turns radius wireGauge_car turns_car radius_car height spacing velocity totalCharge];

totalCharge = DWPTeff(V, wireGauge, turns, radius, wireGauge_car, turns_car, radius_car, height, spacing, velocity, 2);
data(2,:) = [V wireGauge turns radius wireGauge_car turns_car radius_car height spacing velocity totalCharge];

writematrix(data,'data/data.csv')