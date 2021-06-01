% Metas.UncLib.Matlab Tutorial
%
% Quick start guide on how to use Metas.UncLib.Matlab. Provides simple
% examples on how to deal with uncertain objects and points out potential 
% pitfalls.
% It is expected that the user has basic knowledge of MATLAB.
%
% M. Wollensack & M. Zeier, METAS
%

clear all;
echo on;

%% --- LinProp ---
% define function handle for LinProp
unc = @LinProp;

%% Create uncertain numbers
x1 = unc(1,0.1)
x2 = unc(2,0.2)

%% Create uncertain arrays
v1 = unc([1 2],diag([0.1^2 0.2^2]))
% or
v2 = [x1 x2]

%% Create uncertain matrices
% the order in the covariance matrix is based on a column-wise 
% vectorization of the matrix   
m1 = unc([1 2; 3 4],diag([0.1^2 0.3^2 0.2^2 0.4^2]))
% or
m2 = [x1 x2; 3.*x1 2.*x2]

%% Careful: The matrices look the same but they are different 
% m1 and m2 have the same values but different dependencies.
% The elements of m1 are uncorrelated base inputs whereas some of the 
% elements of m2 were constructed from base inputs x1 and x2 and
% are thus correlated. The difference is obvious if one looks at the 
% uncertainties of the inverted matrices

m3 = inv(m1)
m4 = inv(m2)

%% A note on matrices
% direct assignment of a subpart of a matrix does not work

% v3(2,:) = [x1 x2]; %DOES NOT WORK! RETURNS ERROR

% instead preallocate matrix first 
v3 = unc(zeros(2))
% and then assign
v3(2,:) = [x1 x2]

% Preallocation of arrays and matrices is generally recommended because it
% improves performance.

% But careful: Preallocation of an array or matrix as doubles and then 
% filling it with uncertain elements will not work

% v3 = [];
% v3(1) = x1; %DOES NOT WORK! RETURNS DOUBLE

%% Do some calculations with uncertain numbers
x3 = x1 + x2
x4 = x1.*x2
x5 = sqrt(x3+x4)
x6 = x3./x4
x7 = x5.*sin(x6)

%% Return expectation value of result
get_value(x7)
% for LinProp the same as function value
get_fcn_value(x7)

%% Calculate uncertainty of result
get_stdunc(x7)

%% Calculate expanded uncertainty of result
% returns lower and upper limit
get_coverage_interval(x7,0.95)

%% Calculate covariance or correlation of result
get_covariance([x6 x7])
get_correlation([x6 x7])

%% Calculate uncertainty contribution of base inputs
% these are the sensitivity coefficients multiplied with the standard
% uncertainties of the base inputs x1 and x2
get_jacobi(x7)

%% Calculate sensitivities w.r.t. intermediate results
% x3 and x4 are intermediate results
j1 = get_jacobi2(x7,[x3 x4])

% NOTE: The calculation of sensitivities w.r.t. intermediate results
% requires special care and can't be done blindly. There are two
% requirements:
% 1. The vector g = [x3 x4] must be complete in the sense that x7 can be
% calculated as a composite function x7 = f(g(x1,x2))
% e.g. get_jacobi2(x7,x3) would return a wrong result
% 2. The elements of g must be linear independent
% e.g. get_jacobi2(x7,[x3 x4]) with x4 = 2.* x3 would crash
% Generally it is not difficult to satisfy these requirement but it is 
% the responsibility of the user to do so. You will not be warned and 
% the returned result might be wrong, but it can be
% checked by summing up uncertainty contributions and comparing with
% combined uncertainty (see next cell). 

%% Calculate uncertainty contribution of intermediate result
u1 = get_unc_component(x7,[x3 x4])
% this is the same as
j1.*get_stdunc([x3 x4])

% NOTE: the squared sum of these uncertainties is not equal to the 
% combined uncertainty of x7, because x3 and x4 are correlated
sqrt(sum(u1.^2)) % not the same as u(x7)

% but by taking correlation into account it is the same as u(x7)
sqrt(j1*get_covariance([x3 x4])*j1') % same as u(x7)

%% Store an uncertainty object
% stores uncertain objects in xml format maintaining all information
xml_file([x6 x7],'test.xml')

%% Reload a stored uncertainty object
% all information is recovered. 
a1 = unc('test.xml','xml_file');
x6n = a1(1)
x7n = a1(2)

% correlation between x6n and x7n is the same as between x7 and x6
get_correlation([x6n x7n])

% also correlations with respect to any other quantities are the same 
get_correlation([x6n x7n x6 x7 x3])

% NOTE: x6n and x7n are in every aspect identical to x6 and x7. 
% Recovery of full information works even if you shut down Matlab 
% and then restart, as long as the uncertainty objects were stored before 
% shutdown and then reloaded after restart again.
% All the necessary information is stored in the files.

%% How to bridge non-analytical parts in the measurement equation
% If the measurement equation consists of a non-analytical 
% part, as e.g. a numerical method like nonlinear least squares, which 
% returns an output p with value p0 and if the sensitivities j1 j2 of 
% p with respect to the inputs x1 and x2 can be determined somehow, 
% it is possible to define p as an uncertainty number with
% p = unc(p0, [x1 x2], [j1 j2], 'system')

% e.g.
p0 = 5.0
j = [2.1 4.5]
p = unc(p0, [x1 x2], j, 'system')

%% --- DistProp ---
% define function handle for DistProp
unc = @DistProp;

%% Initialize DistProp level
% MaxLevel specifies the maximum order of the taylor expansion
DistPropGlobalMaxLevel(3);

% NOTE: There is no upper limit for MaxLevel but be aware that memory 
% consumption and computational burden increase exponentially with the
% value of MaxLevel.
% MaxLevel=3 is a good value which already provides reasonable
% results for many cases at which LinProp fails.
% 
%% Creation of DistProp objects
% the same as for LinProp
% In the future it will be possible to specify distributions other than
% gaussians

y1 = unc(0,0.1)

%% Calculate with DistProp objects
% LinProp would return an uncertainty of 0
y2 = y1.^2
y3 = cos(y1)

%% Expectation and function value are not the same
get_value(y2)
get_fcn_value(y2)
get_value(y3)
get_fcn_value(y3)

%% Uncertainties of DistProp objects
get_stdunc(y2)
get_stdunc(y3)
get_coverage_interval(y2,0.95)
get_coverage_interval(y3,0.95)

% NOTE: This is based on the assumption that the result distribution is 
% gaussian. It is possible to obtain central moments of the result 
% distribution using get_moment(). This information is not yet used to 
% calculate more realistic coverage intervals.

%%
echo off;

