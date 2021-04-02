clear all, close all, clc

%% FORMULA variables

% Transmitting coil
V = [10 11]; % Voltage of solar panels [V]
wireGauge = [8]; % Wire gauge []
turns = [10]; % Number of turns in the coil []
radius = [10]; % Radius of transmitting coil [m]

% Receiving coil
wireGauge_car = [12]; % Wire gauge []
turns_car = [10]; % Number of turns in the coil []
radius_car = [10]; % Radius of transmitting coil [m]

% Coil orientation
height = [2 30]; % Height of the receiving coil above the transmitting coil [m]
spacing = [3 * radius]; % Distance betewen centers of coils in road [m]
velocity = [1]; % Velocity of car [m/s]

%% Cost Analysis Variables

%% MAIN

variables = 10;
totalScenarios = length(V) * length(wireGauge) * length(turns) * length(radius) * length(wireGauge_car) * length(turns_car) * length(radius_car) * length(height) * length(spacing) * length(velocity)
[A B C D E F G H I J] = ndgrid(V, wireGauge, turns, radius, wireGauge_car, turns_car, radius_car, height, spacing, velocity);

A = reshape(A, [1, totalScenarios]);
B = reshape(B, [1, totalScenarios]);
C = reshape(C, [1, totalScenarios]);
D = reshape(D, [1, totalScenarios]);
E = reshape(E, [1, totalScenarios]);
F = reshape(F, [1, totalScenarios]);
G = reshape(G, [1, totalScenarios]);
H = reshape(H, [1, totalScenarios]);
I = reshape(I, [1, totalScenarios]);
J = reshape(J, [1, totalScenarios]);

data = zeros(totalScenarios, variables + 1); % Initialize data matrix
format shortG % print data with the desired detail

parpool(2);
tic
parfor i = 1 : length(A)
    str = sprintf('%d, %d, %d, %d, %d, %d, %d, %d, %d, %d', A(i), B(i), C(i), D(i), E(i), F(i), G(i), H(i), I(i), J(i));
    disp(str);
    scenarioID = i;
    tic
    totalCharge = DWPTeff(A(i), B(i), C(i), D(i), E(i), F(i), G(i), H(i), I(i), J(i), scenarioID);
    toc
    data(i,:) = [A(i) B(i) C(i) D(i) E(i) F(i) G(i) H(i) I(i) J(i) totalCharge];
end
toc

writematrix(data,'data/data.csv')