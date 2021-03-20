%AUTHOR:    Katherine Kemp (katherine.e.kemp@gmail.com)

function data = main_cost(outputFolder)  
    outputFolder = 'output_cost';
    warning('off');

    %% FORMULA variables
    
    % Transmitting coil
    V = [100 250 450 600]; % Voltage of solar panels [V]
    wireGauge = [6 8]; % Wire gauge []
    turns = [410 460 510 560]; % Number of turns in the coil []
    radius = [1 2 3 4]; % Radius of transmitting coil [m]

    % Receiving coil
    wireGauge_car = [8 12]; % Wire gauge []
    turns_car = [230 260 290 320]; % Number of turns in the coil []
    radius_car = [1.5 1.7 1.9 2.1]; %  Radius of receiving coil [m]

    % Coil orientation
    height = [.15 .3]; % Height of the receiving coil above the transmitting coil [m]
    spacing = [0 .18 .36]; % Distance betewen closest edges of coils in road [m]
    velocity = [30]; % Velocity of car [m/s]

    %% MAIN

    totalScenarios = length(V) * length(wireGauge) * length(turns) * length(radius) * length(radius_car) * length(height) * length(spacing);
    [A, B, C, D, E, F, G] = ndgrid(V, wireGauge, turns, radius, radius_car, height, spacing);

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
    fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', 'scenarioID', 'V', 'wireGauge', 'turns', 'radius', 'wireGauge_car', 'turns_car', 'radius_car', 'height', 'spacing', 'velocity', 'totalCharge', 'cost', 'cost/charge');
    fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', '[]', '[V]', '[]', '[]', '[m]', '[]', '[]', '[m]', '[m]', '[m]', '[m/s]', '[C]', '[$]', '[C/$]');
    fclose(fid);
    
    data =  []; % Initialize data matrix
    format shortG % print data with the desired detail
    mod_val = round(totalScenarios/100);
    parpool();

    parfor i = 1 : totalScenarios
        tic
        fid = fopen(name, 'a');
        returnData = get_cost(A(i), B(i), C(i), D(i), wireGauge_car, turns_car, E(i), F(i), G(i), velocity, i, outputFolder);
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

    delete(gcp('nocreate'))
    fclose('all');
end