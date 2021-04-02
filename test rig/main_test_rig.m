%AUTHOR:    Katherine Kemp (katherine.e.kemp@gmail.com)

function data = main_test_rig(outputFolder)  
    outputFolder = 'output';
    warning('off');

    %% Test rig measurements
    
    OD = 44; % Coil outer diameter [mm]
    ID = 20.5; % Coil inner diameter [mm]
    Deff = (OD^3 - ID^3) / (OD^2 - ID^2) / 3; % Coil effective diameter [mm]
    
    radius_arm = 0.5334; % [m]
    circumference = 2 * pi * radius_arm; % length of outside of track [m]

    coil_count = [16 32 48]; % Number of coils on the track []
   
    %% FORMULA variables
    
    % Transmitting coil
    I = [1 3 5 7 9 11 13]; % Current of source [A]
    R = [19] / 1000; % Coil resistance [Ohms]
    turns = [10]; % Number of turns in the coil []
    radius = [Deff/2] / 1000; % Radius of transmitting coil [m]

    % Receiving coil
    R_car = R; %Coil resistance [Ohms]
    turns_car = turns; % Number of turnmain_test_rigs in the coil []
    radius_car = radius; %  Radius of receiving coil [m]

    % Coil orientation
    height = [0 .05 .1 .15 .2 .25 .3]; % Height of the receiving coil above the transmitting coil [m]
    spacing = (circumference - 2 * radius .* coil_count) ./ coil_count; % Distance betewen closest edges of coils in road [m]
    velocity = [100 * circumference / 60]; % Velocity of car [m/s]

    %% MAIN

    totalScenarios = length(I) * length(R) * length(turns) * length(radius) * length(radius_car) * length(height) * length(spacing);
    [A, B, C, D, E, F, G] = ndgrid(I, R, turns, radius, radius_car, height, spacing);

    A = reshape(A, [1, totalScenarios]);
    B = reshape(B, [1, totalScenarios]);
    C = reshape(C, [1, totalScenarios]);
    D = reshape(D, [1, totalScenarios]);
    E = reshape(E, [1, totalScenarios]);
    F = reshape(F, [1, totalScenarios]);
    G = reshape(G, [1, totalScenarios]);
    
    mkdir(outputFolder);
    name = strcat(outputFolder, '/out.txt');
    fid = fopen(name, 'w');
    fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', 'scenarioID', 'i', 'I', 'R', 'turns', 'radius', 'R_car', 'turns_car', 'radius_car', 'height', 'spacing', 'velocity', 'totalCharge');
    fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', '[]', '[]', '[V]', '[]', '[]', '[m]', '[]', '[]', '[m]', '[m]', '[m]', '[m/s]', '[C]');
    fclose(fid);
    
    data =  []; % Initialize data matrix
    format shortG % print data with the desired detail
    mod_val = round(totalScenarios/100);
    pool = parpool();
    files = 'BSmag_get_B.m';
    addAttachedFiles(pool, files);

    parfor i = 1 : totalScenarios
        tic
        fid = fopen(name, 'a');
        returnData = get_field_test_rig(A(i), B(i), C(i), D(i), R_car, turns_car, E(i), F(i), G(i), velocity, i, outputFolder)
        data = [data; returnData];
        
        for j = 1:length(returnData(:,1))
            for k = 1:length(returnData(j,:))
                fprintf(fid, "%d,", returnData(j, k));
            end
            fprintf(fid, "\n");
        end
        
        if mod(i, mod_val) == 0
            sprintf('%d percent done, i = %d out of %d', round(100*i/totalScenarios), i, totalScenarios)
        end
        toc
    end

    writematrix(data, 'backup.csv')
    delete(gcp('nocreate'))
    fclose('all');
end