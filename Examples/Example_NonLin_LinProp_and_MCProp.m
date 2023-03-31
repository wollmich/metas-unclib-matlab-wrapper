% Example Non-Linear using LinProp and MCProp
% Michael Wollensack METAS - 03.02.2023 - 21.02.2023

clear all;
close all;

disp(sprintf('\nExample Non-Linear\n'))

%% Definition of the inputs using LinProp
x = LinProp(0.0, 1.0, 'x');

%% Compute the outputs using LinProp
y = x.*x;
z = 10.*y;

outputs = [y z]
cv = get_covariance(outputs)

%% Convert LinProp to MCProp
xMC = LinProp2MCProp(x);

%% Compute the outputs using MCProp
yMC = xMC.*xMC;
zMC = 10.*yMC;

outputsMC = [yMC zMC]
cvMC = get_covariance(outputsMC)

%% Convert MCProp to LinProp
outputs2 = MCProp2LinProp(outputsMC, xMC, x)
cv2 = get_covariance(outputs2)

%% Test
test = cv2 - cvMC
