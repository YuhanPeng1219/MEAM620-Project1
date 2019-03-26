close all;
clear all;
clc;
addpath(genpath('./'));

%% Plan path
disp('Planning ...');
map = load_map('map1.txt', 0.2, 0.2, 0.2);
start = {[0.0, -4.9, 0.5]};
stop = {[8.0, 18.0, 3.0]};
%  start = {[2.0, 0.9, 1.2]};
%  stop = {[20.0, 5.0, 6.0]};
%  start = {[0.0, -5.0, 1.2]};
%  stop = {[10.0, 25.0, 5.0]};
%  start = {[5.0, -14, -2]};
%  stop = {[1, 17.0, 2.5]};
nquad = length(start);
for qn = 1:nquad
    path{qn} = dijkstra(map, start{qn}, stop{qn}, true);
end
if nquad == 1
    plot_path(map, path{1});
else
    % you could modify your plot_path to handle cell input for multiple robots
end

%% Generate trajectory
disp('Generating Trajectory ...');
trajectory_generator1([], [], map, path{1});

%% Run trajectory
trajectory = test_trajectory(start, stop, map, path, true); % with visualization
