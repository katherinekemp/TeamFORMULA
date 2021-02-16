%%
%AUTHOR:    20150407, L. Queval (loic.queval@gmail.com)
%COPYRIGHT: 2015, Loic Quéval, BSD License (http://opensource.org/licenses/BSD-3-Clause).

function totalCharge = DWPTz(I, turns, radius, turns_car, radius_car, height, spacing, velocity, scenarioID) 
    %% FORMULA Constants

    maxDistance = 160; % distance the car will travel [m]
    increment = .5; % Resolution (distance between the points in meshgrid) [m]
    distanceStep = 1; % Distance between timesteps [m]
    
    efficiencyOfRectifier = 1; % We assume the rectifier is perfectly efficient for our simulations []
    muRel = 1; % We assume vacuum permeability in our simulations []
    
    dGamma = 1e9; % filament max discretization step [m]
    filamentStep = 100 % number of points in the filament []
    tightness = 10000; % Increasing tightness will make the coil wrapped more tightly [1/m]
    tightness_car = 10000; % Increasing tightness will make the coil wrapped more tightly [1/m]

    %% Initialize
    tic
    plot_scrsz = get(groot,'ScreenSize');
    distance = 0:distanceStep:maxDistance;
    time = distance./velocity;
    totalFLux = zeros(1, length(distance)); % [C]
    BSmag = BSmag_init(); % Initialize BSmag analysis

    %% Initialize video
    %{
    filename1 = sprintf('scenario%u', scenarioID);
    myVideo1 = VideoWriter(filename1); %open video file
    myVideo1.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
    open(myVideo1)

    filename2 = sprintf('scenario%usideview', scenarioID);
    myVideo2 = VideoWriter(filename2); %open video file
    myVideo2.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
    open(myVideo2)
    %}
    %% Create Filaments
    Gamma = [];

    coilCount = (maxDistance - radius) / spacing;
    sign = 1;
    for x = 0:coilCount
        % Source points (points where there is a current source)
        theta = linspace(0,turns*2*pi,turns*filamentStep);
        Gamma = [cos(theta')*radius + radius + x*spacing,sin(theta')*radius,theta'/tightness]; % x,y,z [m,m,m]
        [BSmag] = BSmag_add_filament(BSmag,Gamma,sign*I,dGamma); %% populate BSmag data structure
        sign = -sign;
    end
    
    %% Place vector field
    % Field points (where we want to calculate the field)
    maxRadius = max(radius,radius_car);
    numberOfSquares = 1 + 2*maxRadius / increment; % Calculate the number of squares on the mesh
    x_M = linspace(distance(i),2*maxRadius + distance(i),numberOfSquares); % x [m]
    y_M = linspace(-maxRadius,maxRadius,numberOfSquares); % y [m]
    z_M = linspace(-max(maxRadius, 2*height),max(maxRadius,2*height),numberOfSquares); % z [m]
    [X_M,Y_M,Z_M]=meshgrid(x_M,y_M,z_M);

    %% Biot-Savart Integration
    [X,Y,Z,BZ] = BSmag_get_B(BSmag,X_M,Y_M,Z_M,muRel); %% This is where the integation happens
    BZ(abs(BZ)<1e-10) = 0; % Make tiny values equal to 0

    %% Simulation
    %for i=1:length(distance)
    for i=10
        %% Calculate current flux and total charge
        zeroMark = .5 + numberOfSquares / 2;
        heightIndex = zeroMark + height / increment;

        n = 0;
        flux = zeros(1,numberOfSquares - zeroMark);
        for m = zeroMark:numberOfSquares
            n = n + 1;
            flux(n) = efficiencyOfRectifier*abs(BZ(zeroMark,m,heightIndex));
        end

        totalFlux(i) = 2*pi*increment*trapz(flux);

        %% Car filament for plots
        theta_car = linspace(0,turns_car*2*pi,turns_car*100);
        Gamma_car = [cos(theta_car')*radius_car + radius_car + velocity*time(i),sin(theta_car')*radius_car,theta_car'/tightness_car + height]; % x,y,z [m,m,m]

        %% FIGURE Visualize surface
        %{
        f1 = figure(1);
        f1.OuterPosition = [ 1.5*plot_scrsz(3) .2*plot_scrsz(4) .6*plot_scrsz(3) .7*plot_scrsz(4)];
        hold on, box on, grid on
        plot3(Gamma(:,1),Gamma(:,2),Gamma(:,3),'.-k') % plot filament
        slice(X,Y,Z,BZ,[],[],[height]), colorbar % plot Bz
        plot3(Gamma_car(:,1),Gamma_car(:,2),Gamma_car(:,3),'.-r') % plot car filament
        xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]'), title ('Bz [T]')
        daspect([1,1,1])
        view(3)
        caxis([-1 1]*1e-6)
        drawnow
        pause(.4)

        frame1 = getframe(gcf); % get frame
        writeVideo(myVideo1, frame1);
        clf
        hold off
        %}
        %% FIGURE Visualize surface (side view)
        %{
        f2 = figure(2);
        f2.OuterPosition = [ 1.5*plot_scrsz(3) .2*plot_scrsz(4) .6*plot_scrsz(3) .7*plot_scrsz(4)];
        hold on, box on, grid on
        plot3(Gamma(:,1),Gamma(:,2),Gamma(:,3),'.-k') % plot filament
        slice(X,Y,Z,BZ,[],[0],[]), colorbar % plot Bz
        plot3(Gamma_car(:,1),Gamma_car(:,2),Gamma_car(:,3),'.-r') % plot car filament
        xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]'), title ('Bz [T]')
        daspect([1,1,1])
        view(0,0)
        caxis([-1 1]*1e-6)
        drawnow
        pause(.4)

        frame2 = getframe(gcf); %get frame
        writeVideo(myVideo2, frame2);
        clf
        hold off
        %}
    end

    %close(myVideo1)
    %close(myVideo2)
    inducedEMF = -turns_car*diff(totalFlux)./diff(time);
    
    totalCharge = sum(charge);
    f3 = figure(3);
        hold on
        plot(time,charge)
        xlabel ('time [s]'), ylabel ('Total Accumulated Charge [C]'), title ('Charge vs. Time')
        saveas(f3,sprintf('chargeVStime%u.jpg',scenarioID))
        hold off
    close all
    toc
end