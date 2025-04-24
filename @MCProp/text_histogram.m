function t = text_histogram(x, varargin)
% TEXT_HISTOGRAM Returns the histogram as text.
%
% text_histogram(unc) returns the histogram of unc as text.
% The input unc must be a MCProp. It must be a scalar.
%
% text_histogram(unc, nbins) behaves as above, but also specifies the
% number of bins.
%
% text_histogram(unc, nbins, nlines) behaves as above, but specifies the
% number of lines.
%
% text_histogram(unc, nbins, nlines, fillchar) behaves as above, but
% specifies the fill character.

% Michael Wollensack METAS - 24.04.2025

x = MCProp(x);
n = get_net_object(x);

nbins = 80;
nlines = 24;
fillchar = 'X';

if nargin > 1
    nbins = int32(varargin{1});
end
if nargin > 2
    nlines = int32(varargin{2});
end
if nargin > 3
    fillchar = System.Char(varargin{3});
end
t = n.TextHistogram(nbins, nlines, fillchar);
t = replace(char(t), char([13 10]), newline);
end
