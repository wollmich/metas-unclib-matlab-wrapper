function t = text_unc_budget(x, varargin)
% TEXT_UNC_BUDGET Returns the uncertainy budget as text.
%
% text_unc_budget(unc) returns the uncertainty budget of unc as text.
% The input unc must be a LinProp. It must be a scalar.
%
% text_unc_budget(unc, showId) behaves as above, but also specifies if the
% ID column is shown.
%
% text_unc_budget(unc, showId, showDetail) behaves as above, but specifies
% if the Input Distribution, Input Value, Input Standard Uncertainty and
% Sensitivity columns are shown.
%
% text_unc_budget(unc, showId, showDetail, showDof) behaves as above, but
% specifies if the Dof column is shown.
%
% text_unc_budget(unc, showId, showDetail, showDof, format) behaves as 
% above, but also specifies the format of numerical values. format must be 
% a char array or string scalar. 
% Recommended values are 'f3', 'f6', 'f9', 'e3', 'e6', and 'e9', but
% <a href="https://learn.microsoft.com/en-us/dotnet/standard/base-types/standard-numeric-format-strings">standard</a> and <a
% href="https://learn.microsoft.com/en-us/dotnet/standard/base-types/custom-numeric-format-strings">custom</a> .NET numeric format strings will also work.
%
% text_unc_budget(unc, showId, showDetail, showDof, format, sep) behaves as
% above, but also specifies the column separator.

% Michael Wollensack METAS - 24.04.2025

x = LinProp(x);
n = get_net_object(x);

showId = false;
showDetail = true;
showDof = false;
format = 'f6';
sep = '  ';

if nargin > 1
    showId = logical(varargin{1});
end
if nargin > 2
    showDetail = logical(varargin{2});
end
if nargin > 3
    showDof = logical(varargin{3});
end
if nargin > 4
    format = System.String(varargin{4});
end
if nargin > 5
    sep = System.String(varargin{5});
end
t = n.TextUncBudget(showId, showDetail, showDof, format, sep);
t = replace(char(t), char([13 10]), newline);
end
