%AUTHOR:    Katherine Kemp (katherine.e.kemp@gmail.com)

function totalCharge = get_flux(turns, wireGauge_car, turns_car, radius_car, height, spacing, velocity, scenarioID, BZ, X_M, Y_M, increment, R_car, meshDistance, heightIndex, numberOfSquaresX, numberOfSquaresY, outputFolder) 
    %% FORMULA Constants
    
    maxDistance = 1600; % distance the car will travel [m] 1600, about 1 mile
    distanceStep = increment; % Distance between timesteps [m] keep equal to increment for best results, or multiply by a positive integer for a more coarse simulation
    efficiencyOfRectifier = 1; % We assume the rectifier is perfectly efficient for our simulations []
    tightness_car = 10000; % Increasing tightness will make the coil wrapped more tightly [1/m] We want them to be wrapped perfectly tightly, and this appears to be close enough.

    %% Car and Graph Calculations
    
    %BZlimits = [min(min(BZ(:,:,heightIndex))) max(max(BZ(:,:,heightIndex)))]; % get [minBZ maxBZ]
    
    distance = 0:distanceStep:meshDistance; % Make distance array the max distance - the car coil (so it doesn't go off the page)
    actualDistance = 0:distanceStep:maxDistance;
    time = actualDistance./velocity; % Corresponding time array
    
    flux = zeros(1, length(actualDistance)); % [C] Initialize flux array
        
    %% Initialize video
    %{
    plot_scrsz = get(groot,'ScreenSize'); % Make graphs the screen size
    
    filename1 = sprintf('data/Scenario%u', scenarioID);
    myVideo1 = VideoWriter(filename1, 'MPEG-4'); % open video file
    myVideo1.FrameRate = 10;  % can adjust this, 5 - 10 works well for me
    open(myVideo1)
    %}
    %% Simulation
    
    for i = 2 * spacing / increment : 4 * spacing / increment - 1
        %% Calculate current flux and total charge
        
        minX = distance(i) * numberOfSquaresX / meshDistance; % index in the mesh grid of the minimum x value of the car coil
        field = zeros(numberOfSquaresY, numberOfSquaresY); % initialize flux array to zeros

        for n = round(minX)+1:round(minX)+numberOfSquaresY % Integration accross x
            yLim = sqrt(radius_car^2 - (X_M(1,n) - radius_car - distance(i))^2); % The limits in each row of meshgrid for the y values (so we ignore values outside of the car coil and set them to zero)
            for m = 1:numberOfSquaresY % Integration across y
                if Y_M(m,1) > -1*yLim && Y_M(m,1) < yLim % If the element is within the y limits
                    field(m, n - round(minX)) = BZ(m,n,heightIndex); % put BZ of that location in the field array, we only care about BZ because it is perpendicular
                end
            end
        end

        flux(i) = trapz(Y_M(:,1)',trapz(X_M(1,round(minX)+1:round(minX)+numberOfSquaresY),field,2)); % Integrate field of x and then over y to get total flux at that location

        %% Car filament for plots
        %{
        theta_car = linspace(0, turns_car*2*pi, turns_car*filamentStep);
        Gamma_car = [cos(theta_car')*radius_car + radius_car + distance(i), sin(theta_car')*radius_car, theta_car'/tightness_car + height]; % x,y,z
        %}
        %% FIGURE Visualize surface
        %{
        f1 = figure(1);
        f1.OuterPosition = [1.5*plot_scrsz(3) .2*plot_scrsz(4) .6*plot_scrsz(3) .7*plot_scrsz(4)];
        hold on, box on, grid on
        xlim([0, meshDistance]);
        for n=1:BSmag.Nfilament % Plot all filaments in road in black
            plot3(BSmag.filament(n).Gamma(:,1),BSmag.filament(n).Gamma(:,2),BSmag.filament(n).Gamma(:,3),'.-k') % plot filament
        end
        slice(X,Y,Z,BZ,[],[],[height]), colorbar % plot Bz at heigh of the car coil
        plot3(Gamma_car(:,1),Gamma_car(:,2),Gamma_car(:,3),'.-r') % plot car filament in red
        xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]'), title ('Bz [T] 3D Road View')
        daspect([1,1,1])
        view(3) % Show 3D view
        caxis(1.05*BZlimits)
        drawnow
        pause(.4)

        frame1 = getframe(gcf); % Get frame
        writeVideo(myVideo1, frame1); % Add frame to the video
        clf
        hold off
        %}
    end
    
    k = 4 * spacing / increment - 1;
    for j = 2 * spacing / increment - 1: -1 : 1
        if flux(j) == 0
            flux(j) = flux(k);
        end
        
        if k == 2 * spacing / increment 
            k = 4 * spacing / increment - 1;
        else
            k = k - 1;
        end
    end
    
    k = 2 * spacing / increment;
    for j = 4 * spacing / increment : length(flux)
        if flux(j) == 0
            flux(j) = flux(k);
        end
        
        if k == 4 * spacing / increment - 1
            k = 2 * spacing / increment;
        else
            k = k + 1;
        end
    end
    
    %close(myVideo1)
    
    dT = diff(time);
    inducedEMF = -turns_car*diff(flux)./dT;
    I_car = inducedEMF./R_car;
    totalCharge = trapz(time(1:length(I_car)),efficiencyOfRectifier*abs(I_car));
    cumulativeCharge = cumtrapz(time(1:length(I_car)),efficiencyOfRectifier*abs(I_car));
    
    %{
    f3 = figure(3);
        hold on
        plot(time(1:length(I_car)),I_car)
        xlabel ('time [s]'), ylabel ('Current [A]'), title ('Current vs. Time')
        saveas(f3,fullfile(outputFolder,sprintf('currentVStime%u.jpg',scenarioID)))
        hold off
        
    f4 = figure(4);
        hold on
        plot(time(1:length(cumulativeCharge)),cumulativeCharge)
        xlabel ('time [s]'), ylabel ('Cumulative Charge [C]'), title ('Cumulative Charge vs. Time')
        saveas(f4,fullfile(outputFolder,sprintf('cumChargeVStime%u.jpg',scenarioID)))
        hold off
    %}
    close all
end