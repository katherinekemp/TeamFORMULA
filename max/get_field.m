%AUTHOR:    Katherine Kemp (katherine.e.kemp@gmail.com)

function data = get_field(V, wireGauge, turns, radius, wireGauge_car, turns_car, radius_car, height, original_spacing, velocity, scenarioID, outputFolder) 
    %% FORMULA Constants
    
    increment = .05; % Resolution (distance between the points in meshgrid) [m] 1/50 pi/4 ~= .78375
    muRel = 1; % We assume vacuum permeability in our simulations []
    
    rho = .0171 / 1000^2; % Resistivity of copper [ohm-m]
    density = 8940 * 1000; % Density of copper [g/m^3] https://copperalliance.org.uk/knowledge-base/education/education-resources/bulk-properties-copper-density-resistivity/#:~:text=Density%20is%20the%20mass%20in,(%20kg%20m%2D3).
    wireGauges = [6 8 12]; % Possible wire guages []
    wireDiameters = [.41148 .32639 .20523] / 100; % Corresponding diameters of wire gauges [m] https://ohmslawcalculator.com/awg-wire-chart
    
    maxDistance = 1600; % distance the car will travel [m] 1600, about 1 mile
    dGamma = 1e9; % filament max discretization step [m] found no need to change from provided examples
    filamentStep = 100 * radius; % number of points in the each turn of the filament [1/m]
    tightness = 10000; % Increasing tightness will make the coil wrapped more tightly [1/m] We want them to be wrapped perfectly tightly, and this appears to be close enough.

    %% Initialize
    spacing = 2 * radius + original_spacing;
    
    d_car = zeros(1,length(wireGauge_car));
    for i = 1:length(wireGauge_car)
        d_car(i) = wireDiameters(wireGauges == wireGauge_car(i)); % Wire diameter corresponding to gauge [m]
    end

    d = zeros(1,length(wireGauge));
    for i = 1:length(wireGauge)
        d(i) = wireDiameters(wireGauges == wireGauge(i)); % Wire diameter corresponding to gauge [m]
    end
    A = pi * d^2 / 4; % Cross sectional area [m^2]
    L = 2 * pi * radius * (turns + 1); % Length of a coil [m], we assume the wire beyonf the coil to be negligible in comparison but add one extra loop to get a closer estimate
    R = rho * L / A; % Resistance of a coil [ohms]
    I = V / R; % Currentin road [A]
    
    number_of_coils = maxDistance / spacing;
    cost = number_of_coils * (.6 * I * V + .95 * density * A * L / 453.592); % Cost of the configuration [$]
    meshDistance = round(2 * radius + 6 * spacing);

    BSmag.Nfilament = 0; % Initialize number of source filament (from BSmag_init)

    %% Create Filaments
    
    Gamma = [];
    coilCount = 7;
    sign = 1;
    for x = 0:coilCount
        theta = linspace(0, turns*2*pi, turns*filamentStep); % Source points (points where there is a current source)
        Gamma = [cos(theta')*radius + radius + x * spacing, sin(theta') * radius, theta'/tightness]; % x,y,z
        [BSmag] = BSmag_add_filament(BSmag,Gamma,sign*I,dGamma); %% populate BSmag data structure
        sign = -sign;
    end
    
    %% Place vector field
    % Field points (where we want to calculate the field)
    
    maxRadius = max(radius,radius_car);
    numberOfSquaresX = round(1 + meshDistance / increment); % Calculate the number of squares on the mesh X
    numberOfSquaresY = round(1 + 2*maxRadius / increment); % Calculate the number of squares on the mesh Y
    numberOfSquaresZ = 5; % Calculate the number of squares on the mesh Z
    x_M = linspace(0, meshDistance, numberOfSquaresX); % x [m]
    y_M = linspace(-maxRadius, maxRadius, numberOfSquaresY); % y [m]
    z_M = linspace(height - 2.5 * increment, height + 2.5 * increment, numberOfSquaresZ); % z [m]
    [X_M,Y_M,Z_M]=meshgrid(x_M,y_M,z_M);
    heightIndex = .5 + numberOfSquaresZ / 2;
        
    %% Biot-Savart Integration
    
    [~,~,~,BZ] = BSmag_get_B(BSmag,X_M,Y_M,Z_M,muRel); % This is where the integation happens
    BZ(abs(BZ)<1e-10) = 0; % Make tiny values equal to 0

    %% Flux Integration
    
    totalScenarios = length(d_car) * length(turns_car) * length(velocity);
    [A, B, C] = ndgrid(d_car, turns_car, velocity);
    [D, ~, ~] = ndgrid(wireGauge_car, turns_car, velocity);
    
    A = reshape(A, [1, totalScenarios]);
    B = reshape(B, [1, totalScenarios]);
    C = reshape(C, [1, totalScenarios]);
    D = reshape(D, [1, totalScenarios]);
    
    data = zeros(totalScenarios, 10 + 5); % Initialize data matrix
    format shortG % print data with the desired detail
    
    for i = 1:length(A)
        totalCharge = get_flux(turns, A(i), B(i), radius_car, height, spacing, C(i), rho, (scenarioID - 1)*length(A) + i, BZ, X_M, Y_M, maxDistance, increment, meshDistance, heightIndex, numberOfSquaresX, numberOfSquaresY, outputFolder);
        data(i,:) = [((scenarioID - 1)*length(A) + i) scenarioID V wireGauge turns radius D(i) B(i) radius_car height original_spacing C(i) totalCharge cost cost/totalCharge];
    end
    
end