% Author: Jihua Zhu
% Date: 2017-10-07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
global Model Data
load dragon                                    % load the Model, Data and initial parameters
tic
[R, t, TData] = BdICP(R0, t0, 150);            % Resgistration
toc
%  show the registration results
figure 
plot3(Model(1,:),Model(2,:),Model(3,:),'.r');
hold on;
plot3(TData(1,:),TData(2,:),TData(3,:),'.g');


