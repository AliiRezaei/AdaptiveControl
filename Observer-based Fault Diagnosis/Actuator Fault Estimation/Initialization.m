clc
clear
close all
set(0, 'defaultlinelinewidth', 2)
set(0, 'defaultTextInterpreter', 'latex')

%% NN Dimensions

n = 2;                                    % system order
inputLayerSize  = n + 1;                  % number of neurons in input  layer
hiddenLayerSize = 2 * inputLayerSize + 1; % number of neurons in hidden layer
outputLayerSize = 2;                      % number of neurons in output layer
nnSize = [n, inputLayerSize, hiddenLayerSize, outputLayerSize];