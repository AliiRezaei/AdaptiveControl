clc
clear
close all
set(0, 'defaultlinelinewidth', 2)
set(0, 'defaultTextInterpreter', 'latex')

%% NN Dimension

n               = 2;                      % order of system 
inputLayerSize  = n + 1;                  % number of neurons in input  layer
hiddenLayerSize = 2 * inputLayerSize + 1; % number of neurons in hidden layer
% hiddenLayerSize = 2;                      % number of neurons in hidden layer
outputLayerSize = 2;                      % number of neurons in output layer




