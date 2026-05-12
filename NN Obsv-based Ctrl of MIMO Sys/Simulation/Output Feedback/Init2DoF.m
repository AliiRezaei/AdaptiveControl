clc
clear 
close all

%% NN Dimension

n = 4; % number of system states
p = 2; % number of system inputs/outputs

inputLayerSizeObsv  = 2 * p + n;
hiddenLayerSizeObsv = (2 * inputLayerSizeObsv + 1) * 1;
outputLayerSizeObsv = p;

inputLayerSizeCtrl  = n + p;
hiddenLayerSizeCtrl = (2 * inputLayerSizeCtrl + 1) * 1;
outputLayerSizeCtrl = p;

%% Initial Conditions

% observer network weights init conds
Vobsv0 = 0.1 * randn(hiddenLayerSizeObsv, inputLayerSizeObsv);
Wobsv0 = 0.1 * randn(outputLayerSizeObsv, hiddenLayerSizeObsv);

% controller network weights init conds
Vctrl0 = 0.1 * randn(hiddenLayerSizeCtrl, inputLayerSizeCtrl);
Wctrl0 = 0.1 * randn(outputLayerSizeCtrl, hiddenLayerSizeCtrl);

% system init cond
x0 = [0; 0.5; 0.1; 0];

% observer init cond
x_hat0 = x0([1, 3, 2, 4]);

%% Control Saturation

u_max =   12 * [1; 0.8];
u_min = - 12 * [1; 0.8];

