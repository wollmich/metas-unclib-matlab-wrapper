% Example Resistor Cube
% Michael Wollensack METAS - 18.11.2008

clear all;
close all;

%% Definition of input uncertainty objects R01 ... R12 and U

unc = @LinProp;

U   = 1;

R01 = unc(50,0.1);
R02 = unc(50,0.1);
R03 = unc(50,0.1);
R04 = unc(50,0.1);
R05 = unc(50,0.1);
R06 = unc(50,0.1);
R07 = unc(50,0.1);
R08 = unc(50,0.1);
R09 = unc(50,0.1);
R10 = unc(50,0.1);
R11 = unc(50,0.1);
R12 = unc(50,0.1);

%% Kirchhoff's circuit laws --> linear equation system

Ux = [0 0 0 0 0 0 U U U U U U]';

Rx = [-1   0   0   1   0   0   0   0   1   0   0   0 ;...
       0  -1   0   0   1   1   0   0   0   0   0   0 ;...
       0   0  -1   0   0   0   1   1   0   0   0   0 ;...
       0   0   0   1   1   0   0   0   0  -1   0   0 ;...
       0   0   0   0   0   1   1   0   0   0  -1   0 ;...
       0   0   0   0   0   0   0   1   1   0   0  -1 ;...
      R01  0   0  R04  0   0   0   0   0  R10  0   0 ;...
      R01  0   0   0   0   0   0   0  R09  0   0  R12;...
       0  R02  0   0   0  R06  0   0   0   0  R11  0 ;...
       0  R02  0   0  R05  0   0   0   0  R10  0   0 ;...
       0   0  R03  0   0   0   0  R08  0   0   0  R12;...
       0   0  R03  0   0   0  R07  0   0   0  R11  0 ];

%% Solve linear equation system

Ix = Rx\Ux;

%% Compute equivalent resistor of the cube

I = Ix(1) + Ix(2) + Ix(3);
R = U/I;

%% Output

Ix
R
