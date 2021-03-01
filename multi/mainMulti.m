%AUTHOR:    Katherine Kemp (katherine.e.kemp@gmail.com)

function data = mainMulti() 
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
    velocity = [1 2 3]; % Velocity of car [m/s]

    %% Cost Analysis Variables

    %% MAIN

    variables = 10;
    totalScenarios = length(V) * length(wireGauge) * length(turns) * length(radius) * length(radius_car) * length(height) * length(spacing)
    [A B C D E F G] = ndgrid(V, wireGauge, turns, radius, radius_car, height, spacing);

    A = reshape(A, [1, totalScenarios]);
    B = reshape(B, [1, totalScenarios]);
    C = reshape(C, [1, totalScenarios]);
    D = reshape(D, [1, totalScenarios]);
    E = reshape(E, [1, totalScenarios]);
    F = reshape(F, [1, totalScenarios]);
    G = reshape(G, [1, totalScenarios]);

    data =  []; % Initialize data matrix
    format shortG % print data with the desired detail

    parpool(2);

    tic
    parfor i = 1 : totalScenarios
        %str = sprintf('%d, %d, %d, %d, %d, %d, %d', A(i), B(i), C(i), D(i), E(i), F(i), G(i));
        %disp(str);
        tic
        returnData = get_field(A(i), B(i), C(i), D(i), wireGauge_car, turns_car, E(i), F(i), G(i), velocity, i);
        toc
        data = [data; returnData];
        sprintf('%d percent done, i = %d', 100*i/totalScenarios, i)
    end
    toc

    delete(gcp('nocreate'))
    data
    writematrix(data,fullfile('data','data.csv'))
end