clear all, close all, clc


%% Cost Analysis Variables


%% MAIN

totalScenarios = 2;
data = zeros(totalScenarios, 11); % Initialize data matrix
format shortG % print data with the desired detail

%{
for i = 1:totalScenarios
    scenarioID = 2;
    totalCharge = DWPTeff(V, wireGauge, turns, radius, wireGauge_car, turns_car, radius_car, height, spacing, velocity, scenarioID);
    data(i,:) = [V wireGauge turns radius wireGauge_car turns_car radius_car height spacing velocity totalCharge];
end
%}
totalCharge = DWPT(V, wireGauge, turns, radius, wireGauge_car, turns_car, radius_car, height, spacing, velocity, 1);
data(1,:) = [V wireGauge turns radius wireGauge_car turns_car radius_car height spacing velocity totalCharge];

totalCharge = DWPTeff(V, wireGauge, turns, radius, wireGauge_car, turns_car, radius_car, height, spacing, velocity, 2);
data(2,:) = [V wireGauge turns radius wireGauge_car turns_car radius_car height spacing velocity totalCharge];

writematrix(data,'data/data.out')