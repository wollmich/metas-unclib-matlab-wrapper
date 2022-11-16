function unc_budget(x, varargin)
% UNC_BUDGET Opens budget window.
%
% unc_budget(unc) opens a window with the uncertainty budget of unc.
% The input unc must be a LinProp. It can be a scalar, vector, or matrix.
%
% unc_budget(unc, format) behaves as above, but also specifies the format
% of numerical values. format must be a char array or string scalar. 
% Recomended values are 'f3', 'f6', 'f9', 'e3', 'e6', and 'e9', but
% <a href="https://learn.microsoft.com/en-us/dotnet/standard/base-types/standard-numeric-format-strings">standard</a> and <a
% href="https://learn.microsoft.com/en-us/dotnet/standard/base-types/custom-numeric-format-strings">custom</a> .NET numeric format strings will also work.
%
% unc_budget(unc, format, windowName) behaves as above, but also specifies
% the name of the budget window. windowName must be a string scalar or char
% array. The default value is the name of the variable passed as unc.
%
% unc_budget(unc, format, windowName, columnNames) behaves as above, but
% also specifies the names of the elements of unc and thus the column
% headings in the GUI. columnNames must be a cell string or array of
% strings. This option is only available when var is complex or a vector 
% or matrix.
%
% If unc contains a complex values vector, the real and imag part are
% interleaved. In this respect unc_budget behaves differently than
% get_correlation(), get_covariance(), etc.

% Michael Wollensack METAS - 20.09.2021

x = LinProp(x);
n = get_net_object(x);

name = inputname(1);
if isa(n, 'Metas.UncLib.LinProp.UncNumber')
    c = Metas.UncLib.LinProp.Gui.GuiUncBudget();
else
    c = Metas.UncLib.LinProp.Gui.GuiUncListBudget();
end
try
    c.UnknownIdDecoder = Metas.UncLib.LinProp.UnknownIdDecoderDelegate('Metas.Vna.Data.UncIdDefs','GetInfluenceInfo2');
catch err
end
if nargin > 1
    c.Format = System.String(varargin{1});
end
if nargin > 2
    name = varargin{2};
end
if nargin > 3
    temp = varargin{3};
    n2 = length(temp);
    infos = NET.createArray('System.String',n2);
    for i = 1:n2
        infos(i) = [temp{i} char(10)];
    end
    c.Infos = infos;
end
if isa(n, 'Metas.UncLib.LinProp.UncNumber')
    c.Value = n;
else
    temp = Metas.UncLib.LinProp.UncList();
    l = temp.op_Implicit(n);
    c.Values = l;
end
c.Dock = System.Windows.Forms.DockStyle.Fill;
f = System.Windows.Forms.Form();
f.Text = name;
f.Size.Width = 640;
f.Size.Height = 480;
f.Controls.Add(c);
f.Show();
f.Activate();

% Create a hidden figure and add a callback that is executed on close all.
% Close the form window when the hidden figure is closed by close all.
last_current_figure = get(0, 'CurrentFigure');
fig = figure('visible','off');
fig.UserData = 'Invisible figure to catch ''close all'' and trigger closing of the unc_budget window.';
fig.CloseRequestFcn = {@closeRequest, f};
set(0, 'CurrentFigure', last_current_figure);
end

function closeRequest(fig, ~, form)
% Closes both the Windows.Form and the hidden matlab figure
    form.Close();
    delete(fig);
end
