% LinProp Uncertainty Budget
% Michael Wollensack METAS - 08.08.2017

function unc_budget(x, varargin)

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
fig  = figure('visible','off');
fig.UserData = 'Invisible figure to catch ''close all'' and trigger closing of the unc_budget window.';
fig.CloseRequestFcn = {@closeRequest, f};
end

function closeRequest(~, ~, f)
     f.Close();
end
