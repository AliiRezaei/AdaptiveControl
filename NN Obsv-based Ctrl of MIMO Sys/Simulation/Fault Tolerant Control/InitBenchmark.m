clc
clear 
close all

%% NN Dimension

n   = 3; % number of system states
rho = 2; % relative degree
p   = 1; % number of system inputs/outputs

inputLayerSizeObsv  = 2 * p + rho;
hiddenLayerSizeObsv = 5;
outputLayerSizeObsv = p;

inputLayerSizeCtrl  = rho + p;
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
zeta0 = [0.5; 0; 0.5];

% observer init cond
x_hat0 = [-0.2; 0];

%% Control Saturation

u_max =   2;
u_min = - 2;

q  = 50;

Ts = 1e-3;