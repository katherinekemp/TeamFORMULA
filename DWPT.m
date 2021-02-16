%%
%AUTHOR:    20150407, L. Queval (loic.queval@gmail.com)
%COPYRIGHT: 2015, Loic Quéval, BSD License (http://opensource.org/licenses/BSD-3-Clause).

function totalCharge = DWPTz(I, turns, radius, turns_car, radius_car, height, spacing, velocity, rho_car, L_car, A_car, scenarioID) 
    %% FORMULA Constants

    maxDistance = 160; % distance the car will travel [m] 1600, about 1 mile
    increment = .5; % Resolution (distance between the points in meshgrid) [m] 1/50 pi/4 ~= .78375
    distanceStep = increment; % Distance between timesteps [m] keep equal to increment for best results, or multiply by a positive integer for a more coarse simulation
    
    efficiencyOfRectifier = 1; % We assume the rectifier is perfectly efficient for our simulations []
    muRel = 1; % We assume vacuum permeability in our simulations []
    
    dGamma = 1e9; % filament max discretization step [m] found no need to change from provided examples
    filamentStep = 100; % number of points in the filament [] ???????????
    tightness = 10000; % Increasing tightness will make the coil wrapped more tightly [1/m] a very small number
    tightness_car = 10000; % Increasing tightness will make the coil wrapped more tightly [1/m] a very small number

    %% Initialize
    tic
    plot_scrsz = get(groot,'ScreenSize');
    distance = 0:distanceStep:maxDistance;
    time = distance./velocity;
    totalFlux = zeros(1, length(distance)); % [C]
    BSmag.Nfilament = 0; %Number of source filament

    %% Initialize video
    
    filename1 = sprintf('data/scenario%u', scenarioID);
    myVideo1 = VideoWriter(filename1, 'MPEG-4'); %open video file
    myVideo1.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
    open(myVideo1)

    %{
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
    numberOfSquaresX = 1 + maxDistance / increment; % Calculate the number of squares on the mesh
    numberOfSquaresY = 1 + 2*maxRadius / increment; % Calculate the number of squares on the mesh
    numberOfSquaresZ = 1 + height / increment; % Calculate the number of squares on the mesh
    x_M = linspace(0,maxDistance + 2*maxRadius,numberOfSquaresX + 2*maxRadius / increment + 1); % x [m]
    y_M = linspace(-maxRadius,maxRadius,numberOfSquaresY); % y [m]
    z_M = linspace(.5*height,1.5*height,numberOfSquaresZ); % z [m]
    [X_M,Y_M,Z_M]=meshgrid(x_M,y_M,z_M);
    heightIndex = .5 + numberOfSquaresZ / 2;
    
    %% Biot-Savart Integration
    [X,Y,Z,BZ] = BSmag_get_B(BSmag,X_M,Y_M,Z_M,muRel); %% This is where the integation happens
    BZ(abs(BZ)<1e-10) = 0; % Make tiny values equal to 0
    BZlimits = [min(min(BZ(:,:,heightIndex))) max(max(BZ(:,:,heightIndex)))];
    
    %% Simulation
    %for i=1:length(distance)
    for i=10
        %% Calculate current flux and total charge
        minX = distance(i) * numberOfSquaresX / maxDistance;
        flux = zeros(numberOfSquaresY, numberOfSquaresY);

        for n = round(minX)+1:round(minX)+numberOfSquaresY
            yLim = sqrt(radius_car^2 - (X_M(1,n) - radius_car - distance(i))^2);
            for m = 1:numberOfSquaresY
                if Y_M(m,1) > -1*yLim & Y_M(m,1) < yLim
                    flux(m, n - round(minX)) = BZ(m,n,heightIndex);
                end
            end
        end

        totalFlux(i) = trapz(Y_M(:,1)',trapz(X_M(1,round(minX)+1:round(minX)+numberOfSquaresY),flux,2));

        %% Car filament for plots
        theta_car = linspace(0,turns_car*2*pi,turns_car*filamentStep);
        Gamma_car = [cos(theta_car')*radius_car + radius_car + distance(i),sin(theta_car')*radius_car,theta_car'/tightness_car + height]; % x,y,z [m,m,m]

        %% FIGURE Visualize surface
        
        f1 = figure(1);
        f1.OuterPosition = [1.5*plot_scrsz(3) .2*plot_scrsz(4) .6*plot_scrsz(3) .7*plot_scrsz(4)];
        hold on, box on, grid on
        xlim([0, maxDistance])
        for n=1:BSmag.Nfilament
            plot3(BSmag.filament(n).Gamma(:,1),BSmag.filament(n).Gamma(:,2),BSmag.filament(n).Gamma(:,3),'.-k') % plot filament
        end
        %sliceX = linspace(distance(i) + radius_car,distance(i) + 3*radius_car,numberOfSquaresY);
        slice(X,Y,Z,BZ,[],[],[height]), colorbar % plot Bz
        plot3(Gamma_car(:,1),Gamma_car(:,2),Gamma_car(:,3),'.-r') % plot car filament
        xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]'), title ('Bz [T]')
        daspect([1,1,1])
        view(3)
        caxis(1.05*BZlimits)
        drawnow
        pause(.4)

        frame1 = getframe(gcf); % get frame
        writeVideo(myVideo1, frame1);
        clf
        hold off
        
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

    close(myVideo1)
    %close(myVideo2)
    R_car = rho_car * L_car / A_car;
    dT = diff(time);
    inducedEMF = -turns_car*diff(totalFlux)./dT;
    I_car = inducedEMF./R_car;
    totalCharge = trapz(time(1:length(I_car)),efficiencyOfRectifier*abs(I_car));
    cumCharge = cumtrapz(time(1:length(I_car)),efficiencyOfRectifier*abs(I_car));
    
    f3 = figure(3);
        hold on
        plot(time(1:length(I_car)),I_car)
        xlabel ('time [s]'), ylabel ('Current [A]'), title ('Current vs. Time')
        saveas(f3,sprintf('data/currentVStime%u.jpg',scenarioID))
        hold off
        
    f4 = figure(4);
        hold on
        plot(time(1:length(cumCharge)),cumCharge)
        xlabel ('time [s]'), ylabel ('Cumulative Charge [C]'), title ('Cumulative Charge vs. Time')
        saveas(f4,sprintf('data/cumChargeVStime%u.jpg',scenarioID))
        hold off
        
    close all
    toc
end