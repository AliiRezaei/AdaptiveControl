clc
clear 
close all

%% NN Dimension

n = 6; % number of system states
p = 3; % number of system inputs/outputs

inputLayerSizeObsv  = 2 * p + n;
hiddenLayerSizeObsv = 3;
outputLayerSizeObsv = p;

inputLayerSizeCtrl  = n + p;
hiddenLayerSizeCtrl = 3;
outputLayerSizeCtrl = p;

%% Initial Conditions

% observer network weights init conds
Vobsv0 = 0.1 * randn(hiddenLayerSizeObsv, inputLayerSizeObsv);
Wobsv0 = 0.1 * randn(outputLayerSizeObsv, hiddenLayerSizeObsv);

% controller network weights init conds
Vctrl0 = 0.1 * randn(hiddenLayerSizeCtrl, inputLayerSizeCtrl);
Wctrl0 = 0.1 * randn(outputLayerSizeCtrl, hiddenLayerSizeCtrl);

% system init cond
x0 = zeros(n, 1);

% observer init cond
x_hat0 = zeros(n, 1);

%% Control Saturation

u_max =   25 * [1; 1; 1];
u_min = - 25 * [1; 1; 1];
