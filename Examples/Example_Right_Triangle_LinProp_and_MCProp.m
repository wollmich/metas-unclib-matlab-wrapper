% Example Right Triangle using LinProp and MCProp
% Michael Wollensack METAS - 03.02.2023 - 21.02.2023

clear all;
close all;

disp(sprintf('\nExample Right Triangle\n'))

%% Definition of the inputs using LinProp
a = LinProp(3.0, 0.3, 'a');
b = LinProp(4.0, 0.4, 'b');
inputs = [a b];

%% Compute the outputs using LinProp
c = sqrt(a.*a + b.*b);
p = a + b + c;
A = a.*b./2;

outputs = [c p A]
cv = get_covariance(outputs)

%% Convert LinProp to MCProp
inputsMC = LinProp2MCProp(inputs);

%% Compute the outputs using MCProp
aMC = inputsMC(1);
bMC = inputsMC(2);
cMC = sqrt(aMC.*aMC + bMC.*bMC);
pMC = aMC + bMC + cMC;
AMC = aMC.*bMC./2;

outputsMC = [cMC pMC AMC]
cvMC = get_covariance(outputsMC)

%% Convert MCProp to LinProp
outputs2 = MCProp2LinProp(outputsMC, inputsMC, inputs)
cv2 = get_covariance(outputs2)

%% Test
test = cv2 - cvMC
