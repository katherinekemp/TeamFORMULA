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
clear all, close all, clc

%% FORMULA variables

I = 10; % filament current [A]
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
spacing = 30; % distance betewen centers of coils
velocity = 1; % velocity of car [m/s]


%% Initialize
plot_scrsz = get(groot,'ScreenSize');
time = 0:1:(4*spacing + 2*radius); %% Must start at 0
charge = 0; % [C]

%% Initialize video
filename = 'myVideoFile';
delete filename
myVideo = VideoWriter(filename); %open video file
myVideo.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
open(myVideo)

%% Create Filament
Gamma = [];
for x = 0:4
    % Source points (points where there is a current source)
    theta = linspace(0,turns*2*pi,turns*100);
    Gamma = [Gamma;[cos(theta')*radius + radius + x*spacing,sin(theta')*radius,theta'/tightness]]; % x,y,z [m,m,m]
end

dGamma = 1e9; % filament max discretization step [m]
BSmag = BSmag_init(); % Initialize BSmag analysis
[BSmag] = BSmag_add_filament(BSmag,Gamma,I,dGamma); %% populate BSmag data structure

%% Simulation
for i=1:length(time)
    %% Place vector field
    % Field points (where we want to calculate the field)
    numberOfSquares = 1 + 2*radius / increment; % Calculate the number of squares on the graph
    x_M = linspace(0 + time(i)*velocity,2*radius + time(i)*velocity,numberOfSquares); % x [m]
    y_M = linspace(-radius,radius,numberOfSquares); % y [m]
    z_M = linspace(-radius,radius,numberOfSquares); % z [m]
    [X_M,Y_M,Z_M]=meshgrid(x_M,y_M,z_M);
    
    %BSmag_plot_field_points(BSmag,X_M,Y_M,Z_M); % shows the field points volume %Sets up meshgrid

    %% Biot-Savart Integration
    [BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,X_M,Y_M,Z_M,muRel); %% This is where the integation happens

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
   
    %% FIGURE Visualize surface
    Gamma_car = [cos(theta')*radius + radius + velocity*time(i),sin(theta')*radius,theta'/tightness + height]; % x,y,z [m,m,m]
    
    f = figure(1);
    f.OuterPosition = [ 1.5*plot_scrsz(3) .2*plot_scrsz(4) .6*plot_scrsz(3) .7*plot_scrsz(4)];
    hold on, box on, grid on
    plot3(Gamma(:,1),Gamma(:,2),Gamma(:,3),'.-k') % plot filament
    slice(X,Y,Z,BZ,[],[],[height]), colorbar % plot Bz
    plot3(Gamma_car(:,1),Gamma_car(:,2),Gamma_car(:,3),'.-r') % plot filament
    xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]'), title ('Bz [T]')
    daspect([1,1,1])
    view(3)
    caxis([-1 1]*1e-6)
    drawnow
    pause(.4)
    
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
    clf
    hold off
    
end

close(myVideo)

totalCharge = sum(charge)
plot(time,charge)