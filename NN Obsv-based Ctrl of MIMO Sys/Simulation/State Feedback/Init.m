clc
clear 
close all

%% NN Dimension

n = 4; % number of system states
p = 2; % number of system inputs/outputs

inputLayerSizeCtrl  = n + p;
hiddenLayerSizeCtrl = (2 * inputLayerSizeCtrl + 1) * 1;
outputLayerSizeCtrl = p;

%% Initial Conditions

% controller network weights init conds
Vctrl0 = 0.1 * randn(hiddenLayerSizeCtrl, inputLayerSizeCtrl);
Wctrl0 = 0.1 * randn(outputLayerSizeCtrl, hiddenLayerSizeCtrl);

% system init cond
x0 = [0; 0.5; 0.1; 0];

%% Control Saturation

u_max =   12 * [1; 0.8];
u_min = - 12 * [1; 0.8];

