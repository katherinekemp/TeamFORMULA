%% Current issues:
%{
-Need value of reissitivity
-We can calculate the flux through the coil but we do not know how to
calculate the actual energy in joules that is stored
-Need to turn slice into a circle
%}

%%
%AUTHOR:    20150407, L. Queval (loic.queval@gmail.com)
%COPYRIGHT: 2015, Loic Quéval, BSD License (http://opensource.org/licenses/BSD-3-Clause).

%% Initialize
clear all, close all, clc
%addpath('~/Desktop/Team Formula/Magnetic Fields/BSmag Core/','-end')
plot_scrsz = get(groot,'ScreenSize');
time = 0:pi/10:(2*pi); %% Must start at 0
charge = 0; % [C]

%% Simulation
for i=1:length(time)
    %% FORMULA variables
    I = 10 * sin(time(i)); % filament current [A]
    increment = .5; % Resolution (distance between the points) [m]
    turns = 10; % Numer of turn in the coil []
    radius = 10; % Radius of transmitting coil [mm?]
    tightness = 10000; % Increasing tightness will make the coil wrapped more tightly [1/mm?]
    omega = 1; % Angular frequency of transmitting current [rad/s]
    amplitude = 10; % Amplitude of transmitting current [A]
    muRel = 1; % Relative permeability of material [N/A^2?]
    height = 2; % Height of the receiving coil above the transmitting coil MUST BE ON THE GRID
    efficiencyOfRectifier = 1; % BETWEEN 0 and 1. Approximates the efficiency of the full wave rectifier
    rho = 10; % resistivity of the wire
    distance = 5; % distance betewen centers of coils

    %% Create Filament
    for x = 0:4
        % Source points (points where there is a current source)
        theta = linspace(0,turns*2*pi,turns*100);
        Gamma = [cos(theta')*radius + x*distance,sin(theta')*radius,theta'/tightness]; % x,y,z [m,m,m]

        dGamma = 1e9; % filament max discretization step [m]
        BSmag = BSmag_init(); % Initialize BSmag analysis
        [BSmag(x)] = BSmag_add_filament(BSmag,Gamma,I,dGamma); %% populate BSmag data structure
    end

    %% Place vector field
    % Field points (where we want to calculate the field)
    numberOfSquares = 1 + 2*radius / increment; % Calculate the number of squares on the graph
    x_M = linspace(-radius,radius,numberOfSquares); % x [m]
    y_M = linspace(-radius,radius,numberOfSquares); % y [m]
    z_M = linspace(-radius,radius,numberOfSquares); % z [m]
    [X_M,Y_M,Z_M]=meshgrid(x_M,y_M,z_M);
    
    %BSmag_plot_field_points(BSmag,X_M,Y_M,Z_M); % shows the field points volume %Sets up meshgrid

    %% Biot-Savart Integration
    [BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,X_M,Y_M,Z_M, muRel); %% This is where the integation happens

    BX( abs(BX)<1e-10 ) = 0 ;
    BY( abs(BY)<1e-10 ) = 0 ;
    BZ( abs(BZ)<1e-10 ) = 0 ;
    
    %% Calculate current flux and total charge
    zeroMark = .5 + numberOfSquares / 2;
    heightIndex = zeroMark + height / increment;
    
    n = 0;
    flux = zeros(1,numberOfSquares - zeroMark);
    for m = zeroMark:numberOfSquares
        n = n + 1;
        flux(n) = efficiencyOfRectifier*abs(BZ(zeroMark,m,heightIndex));
    end
    
    charge(i) = 2*pi*increment*trapz(flux);
    
    %% FIGURE Vector Field Plot
    %Plot B/|B|
    
    %{ figure(1)

    max(max(max(abs(BX))));
    max(max(max(abs(BY))));
    max(max(max(abs(BZ))));
    
    f.OuterPosition = [ 1.5*plot_scrsz(3) .2*plot_scrsz(4) .6*plot_scrsz(3) .7*plot_scrsz(4)];
    normB=sqrt(BX.^2+BY.^2+BZ.^2);
    quiver3(X,Y,Z,BX./normB,BY./normB,BZ./normB,'b');
    quiver3(X,Y,Z,1e6*BX,1e6*BY,1e6*BZ,'b');
    drawnow; pause(.5);
    xlim([-10 10]); ylim([-10 10]); zlim([-10 10]); axis tight;
    
   
    %% FIGURE Visualize surface
    
    f = figure(2);
    f.OuterPosition = [ 1.5*plot_scrsz(3) .2*plot_scrsz(4) .6*plot_scrsz(3) .7*plot_scrsz(4)];
    hold on, box on, grid on
    	plot3(Gamma(:,1),Gamma(:,2),Gamma(:,3),'.-r') % plot filament
    	slice(X,Y,Z,BZ,[],[],[height]), colorbar % plot Bz
        %%slice(X,Y,Z,BZ,[0],[],[-radius/2,0,radius/2]), colorbar % plot Bz
    xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]'), title ('Bz [T]')
    view(3), axis equal, axis tight
    caxis([-1 1]*1e-6)
    drawnow
    pause(.4)
    
    
    %% FIGURE Flux tubes
    %{
    f = figure(3);
     f.OuterPosition = [ 1.5*plot_scrsz(3) .2*plot_scrsz(4) .6*plot_scrsz(3) .7*plot_scrsz(4)];
             clf
     figure(3), 
     hold on, box on, grid on
         plot3(Gamma(:,1),Gamma(:,2),Gamma(:,3),'.-r') % plot filament
         [X0,Y0,Z0] = ndgrid(-radius:2:radius,-radius:2:radius,[-radius/2 0 radius/2]); % define tubes starting point        
         htubes = streamtube(stream3(X,Y,Z,BX,BY,BZ,X0,Y0,Z0), [.2 10]);
         streamtube(stream3(X,Y,Z,-BX,-BY,-BZ,X0,Y0,Z0), [.2 10]);
     xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]'), title ('Some flux tubes')
     view(3), axis equal, axis tight
     set(htubes,'EdgeColor','none','FaceColor','c') % change tube color
     camlight left % change tube light
     xlim([-10 10])
     ylim([-10 10])
     zlim([-10 10])
    drawnow
     pause(.4)
    %}
 
    %%
    
end

totalCharge = sum(charge)
plot(time,charge)