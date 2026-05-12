clc
clear
close all
set(0, 'defaultlinelinewidth', 2)
set(0, 'defaultTextInterpreter', 'latex')

%% NN Dimension

n = 2;
p = 1;

inputLayerSizeObsv  = n + p;
% hiddenLayerSizeObsv = 2 * inputLayerSizeObsv + 1;
hiddenLayerSizeObsv = 3;
outputLayerSizeObsv = n;

inputLayerSizeFEstim  = n + p;
% hiddenLayerSizeFEstim = 2 * inputLayerSizeFEstim + 1;
hiddenLayerSizeFEstim = 3;
outputLayerSizeFEstim = p;

Vobsv0 = 0.1 * randn(hiddenLayerSizeObsv, inputLayerSizeObsv);
Wobsv0 = 0.1 * randn(outputLayerSizeObsv, hiddenLayerSizeObsv);

Vfestim0 = 0.1 * randn(hiddenLayerSizeFEstim, inputLayerSizeFEstim);
Wfestim0 = 0.1 * randn(outputLayerSizeFEstim, hiddenLayerSizeFEstim);



