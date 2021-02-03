clear all, close all, clc

%% FORMULA variables

% Transmitting coil
I = 10; % filament current [A]
turns = 10; % Numer of turn in the coil []
radius = 10; % Radius of transmitting coil [mm?]
tightness = 10000; % Increasing tightness will make the coil wrapped more tightly [1/mm?] THIS COULD BE CONSTANT?
rho = 10; % resistivity of the wire

% Receiving coil
turns_car = 10; % Numer of turn in the coil []
radius_car = 10; % Radius of transmitting coil [mm?]
tightness_car = 10000; % Increasing tightness will make the coil wrapped more tightly [1/mm?] THIS COULD BE CONSTANT?
rho_car = 10; % resistivity of the wire

% Coil orientation
height = 2; % Height of the receiving coil above the transmitting coil MUST BE ON THE GRID
spacing = 30; % distance betewen centers of coils
velocity = 1; % velocity of car [m/s]

%% FORMULA Constants

distance = 160; % distance the car will travel [m]
increment = .5; % Resolution (distance between the points) [m]
efficiencyOfRectifier = 1; % BETWEEN 0 and 1. Approximates the efficiency of the full wave rectifier. We will assume it is 1.
muRel = 1; % Relative permeability of road material + air combo [N/A^2?]
scenarioID = 1; % scenario #
timeStep = 1; % Distance between timesteps [s]
dGamma = 1e9; % filament max discretization step [m]

%% MAIN

totalScenarios = 2;
data = zeros(totalScenarios, 13);
format shortG

for i = 1:totalScenarios
    totalCharge = DWPTz(I, turns, radius, tightness, rho, turns_car, radius_car, tightness_car, rho_car, height, spacing, velocity, distance, increment, efficiencyOfRectifier, muRel,scenarioID, timeStep, dGamma);
    data(i,:) = [I turns radius tightness rho turns_car radius_car tightness_car rho_car height spacing velocity totalCharge];
end

writematrix(data,'data.csv')