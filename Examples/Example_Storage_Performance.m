% Example Storage Performance
% Michael Wollensack METAS - 05.07.2012

clear all;
close all;

unc = @LinProp;

n = 20;
fp_bin = 'testdata.bin';
fp_xml = 'testdata.xml';

disp(sprintf('\nExample Storage Performance\n'))

%% Create Random Data
data0 = unc(rand(n, n), eye(n.*n));
data = inv(data0);
disp(sprintf('Memory Size of Data : %.0f kBytes\n', memsize(data)./1024))

%% Binary Serialization
disp(sprintf('Binaray serialization - Save file'))
tic; binary_file(data, fp_bin); toc
disp(sprintf('Binaray serialization - Load file'))
tic; bin_data = unc(fp_bin, 'binary_file'); toc
d = dir(fp_bin);
disp(sprintf('Binaray serialization - Filesize : %.0f kBytes\n', d.bytes./1024))
clear d;
delete(fp_bin);

%% Xml Serialization
disp(sprintf('Xml serialization - Save file'))
tic; xml_file(data, fp_xml); toc
disp(sprintf('Xml serialization - Load file'))
tic; xml_data = unc(fp_xml, 'xml_file'); toc
d = dir(fp_xml);
disp(sprintf('Xml serialization - Filesize : %.0f kBytes\n', d.bytes./1024))
clear d;
delete(fp_xml);
