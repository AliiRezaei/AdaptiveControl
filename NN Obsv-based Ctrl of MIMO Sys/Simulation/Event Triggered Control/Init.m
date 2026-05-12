clc
clear 
close all

%% NN Dimension

n = 2; % number of system states
p = 1; % number of system inputs/outputs

inputLayerSizeObsv  = 2 * p + n;
hiddenLayerSizeObsv = 5;
outputLayerSizeObsv = p;

inputLayerSizeCtrl  = n + p;
hiddenLayerSizeCtrl = 20;
outputLayerSizeCtrl = p;

%% Initial Conditions

% observer network weights init conds
Vobsv0 = 0.1 * randn(hiddenLayerSizeObsv, inputLayerSizeObsv);
Wobsv0 = 0.1 * randn(outputLayerSizeObsv, hiddenLayerSizeObsv);

% controller network weights init conds
Vctrl0 = 0.1 * randn(hiddenLayerSizeCtrl, inputLayerSizeCtrl);
Wctrl0 = 0.1 * randn(outputLayerSizeCtrl, hiddenLayerSizeCtrl);

% system init cond
x0 = [0.5; 0];

% observer init cond
x_hat0 = [0.9; 0];

%% Control Saturation

u_max =   5;
u_min = - 5;
Ts    = 1e-3;