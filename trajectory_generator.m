function [ desired_state ] = trajectory_generator(t, qn, map, path)
% TRAJECTORY_GENERATOR: Turn a Dijkstra or A* path into a trajectory
%
% NOTE: This function would be called with variable number of input
% arguments. At first, it will be called with arguments
% trajectory_generator([], [], map, path) and later, in test_trajectory,
% it will be called with only t and qn as arguments, so your code should
% be able to handle that. This can be done by checking the number of
% arguments to the function using the "nargin" variable, check the
% MATLAB documentation for more information.
%
% parameters:
%   map: The map structure returned by the load_map function
%   path: This is the path re turned by your planner (dijkstra function)
%   desired_state: Contains all the information that is passed to the
%                  controller, as in Phase 1
%
% NOTE: It is suggested to use "persistent" variables to store map and path
% during the initialization call of trajectory_generator, e.g.
% persistent map0 path0
% map0 = map;
% path0 = path;

%% ******************determine the input condition*********************
% Do: Set a varaible to represent the input condition;
%     If the input number is 4, the function is first called 
%     during simulation;
%     If the input number is 2, the function take T and QN as input
%     to generate desire states for controller

if nargin < 4 
   actual_trajectory = true;
else
   actual_trajectory = false;
end

%% *********************Initial varaibles here*************************
% Do: persistent time_matrix restore the time of each waypoint
%     persistent path_refine restore the new path
persistent time_matrix;
persistent path_processed;
persistent snap_coeffs;
persistent total_dist;
persistent dist_matrix;
persistent time_segment;
persistent Total_time;

%% **************Process the trajectory during first call****************
% When variable actual_trajectory = false, take MAP and PATH as input
% to refine the path;
% DO: Check if more than three points are on the same line, delete the 
%     intermediate points. To decrease the waypoints in interpolation;
%     Generate the interpolation matrix and compute the coefficient.
%     Minimum_snap Trajectory interpolation is used here;
%
if actual_trajectory == false
   % ***************** Delete intermediate points *************************
   path_refine = refinePath(path);

   % ************************* Path smooth ********************************
    path_smooth = pathSmooth(path_refine, map);
    path_n=path_smooth;
    
%     resxy_persistent = map.res_xyz(1);
%     points_path = path_refine_samedirection_new(path_refine, 0.1);
%     points_path = path_refine_closepts(points_path, 0, 1.42*resxy_persistent);
%     points_path = path_refine_samedirection_new(points_path, 0.2);
%     points_path(end,:) = path(end,:);
%     path_n = points_path;
   
   % **************** Allocate time to each waypoints *********************
   dist_matrix = sqrt(sum(diff(path_n).^2, 2));
   total_dist = sum(dist_matrix);  
   if total_dist < 15
      Total_time = 9;
   elseif total_dist > 15 && total_dist < 32
       Total_time = 11.5;
   elseif total_dist > 32 && total_dist < 36
       Total_time = 23;
   elseif total_dist > 36 && total_dist < 42
       Total_time = 18;
   elseif total_dist > 42
       Total_time = 23;
   end
   time_segment = Total_time*dist_matrix/total_dist;
   time_matrix = [0;cumsum(time_segment)];
   path_processed = path_n;
   
   % **************** Generate interpolation matrixs* *********************
   for i = 1:size(time_matrix)-1
       coeffs = Min_SnapSpline(0, time_matrix(i+1), ...
                         path_processed(i,:), path_processed(i+1,:));
       snap_coeffs = [snap_coeffs; coeffs];
   end
   
else
    %% ************************ Output assign value************************
%     if t < Total_time
%         [pos, vel, acc] = computeTrajectory(time_matrix, t, snap_coeffs);
% 
%         desired_state.pos = pos';
%         desired_state.vel = vel';
%         desired_state.acc = acc';
%         desired_state.yaw = 0;
%         desired_state.yawdot = 0;   
    t_fregment = time_matrix;
    [n_floor,~] = find(t_fregment <= t);
    n_floor_num = n_floor(end);  
    
    n_ceil_num = n_floor_num + 1;

    if n_ceil_num <= size(path_processed, 1)
        [pos, vel, acc] = traject_interpolation(t_fregment(n_floor_num), t_fregment(n_ceil_num), t, path_processed(n_floor_num,:), path_processed(n_ceil_num,:));
        desired_state.pos = pos(:);  %left: 3x1, right: pos 3x1
        desired_state.vel = vel(:);
        desired_state.acc = acc(:);
                desired_state.yaw = 0;
        desired_state.yawdot = 0;  
    else
        desired_state.pos = path_processed(end, :)';
        desired_state.vel = [0,0,0]';
        desired_state.acc = [0,0,0]';
        desired_state.yaw = 0;
        desired_state.yawdot = 0;  
    end
end


end
%% **************Function to delete useless waypoints*****************
function newpath = refinePath(path)
% REFINE PATH
% function newpath = refinePath(path)
% Do: Given a path with many many waypoints on the same line, delete those
%     intermediate points. Keep only the start and end points of each
%     segment.
%
%     Input:
%       path: [n, 3] matrix. Each row represents a waypoint
%     Output:
%       newpath: [m, 3] matrix. Where m is lower or equal to n.
%
% ***********find the direction between each pair of waypoints*********
newpath = path;
diff_path = diff(newpath);
norm_row = sum(diff_path.^2,2);
path_direct = diff_path./ norm_row;
row_delete = [];

% *************** restore the rows should be canceled *****************
for i = 1:size(path_direct) - 1
    if path_direct(i, :) == path_direct(i+1, :)
        row_delete = [row_delete; i+1];
    end
end

% ********************** delete those rows ****************************
newpath(row_delete, :) = [];
end

%% ************** Function to check line block collision***************
function flag = lineBlock_intersect(start, goal, map)
% LINEBLOCK INTERSECT
% Do: Check whether the line intersects with blocks. If check the function
%     will return a flag.
%     Input:
%       map: a structure with the position and size of blocks
%       start: the start point of the line segment
%       end: the end point of the line segment
%     Output:
%       flag: return 1 if check occurs, else return 0
%
collision = [];
margin = 0.2;
newblocks = zeros(size(map.blocks, 1), 6);
obstacle = map.blocks(:, 1:6);

% ************ To make sure free of collision, add margin to blocks********
newblocks(:, 1:3) = obstacle(:, 1:3) - margin;
newblocks(:, 4:6) = obstacle(:, 4:6) + margin;

% ************************* check every blocks*****************************
for i=1:size(map.blocks)
    c_flag = detectcollision(start, goal, newblocks(i,:));
    collision = [collision,c_flag];
end

% ************************* Return value **********************************
if any(collision == 1)
    flag = 1;
else
    flag = 0;
end
end

%% ****************** Function to smooth path *************************
function path_smooth = pathSmooth(path, map)
% PATH SMOOTH
% Do: Try to line up all the waypoints from start to end. If the two way
%     points can be lined without collision with blocks, write down those
%     waypoints. In this way, the size of path will decrease greatly.
%     Input:
%       map: a structure with the position and size of blocks
%       path: a [n, 3] matrix, in which all points are restored
%     Output:
%       path_smooth: a new path with much less waypoints. [m, 3] matrix,
%       where m < n
%
path_smooth = path(1, :);
size_path = size(path, 1);
i = 1;

while ( ~isequal(path_smooth(i, :), path(end, :)))      
% *********************** Line and check **********************************
    for j = 1:size_path
        if lineBlock_intersect(path_smooth(i,:), path(size_path-j+1, :),...
            map) == 0         % line the two points if no collision
            path_smooth = [path_smooth; path(size_path-j+1, :)];
            break;
        end
    end
    
% ************** Check until the new waypoint is goal *********************
    if path(size_path-j+1, :) == path(end, :)
        break;
    end
    i=i+1;
end
end

%% ****************** Function to generate Spline matrix ******************
function snap_coeffs = Min_SnapSpline(t, tf, start, goal)
% MINIMUM SNAP TRAJECTORY
% DO: Generate a Spline Matrix that can be used to minimize the snap of
%     drone. The trajectory takes the whole path as input and add
%     constrains on intermediate points. As a result, the drone can move
%     without stop and go.
%     Input:
%       t_matrix: A [n, 3] matrix restores the time of each waypoint
%       path: [n, 3] matrix restores the position of each waypoint
%       state_SF: a [6, 3] matrix restores the start and final states. The
%                 first three rows are the vel, acc, snap, for start.
%       time: the time now, used to compute the trajectory in realtime
%     Output:
%       pos: the position at each time [1, 3]
%       vel: the velocity at each time [1, 3]
%       acc: the acceleration at each time [1, 3]
%
% ****************** Write sanp_matrix and state matrix *******************
% snap_matrix = [t^7, t^6, t^5, t^4, t^3, t^2, t, 1;
%               7*t^6, 6*t^5, 5*t^4, 4*t^3, 3*t^2, 2*t, 1, 0;
%               42*t^5, 30*t^4, 20*t^3, 12*t^2, 6*t, 2, 0, 0;
%               210*t^4, 120*t^3, 60*t^2, 24*t, 6, 0, 0, 0;
%               tf^7, tf^6, tf^5, tf^4, tf^3, tf^2, tf, 1;
%               7*tf^6, 6*tf^5, 5*tf^4, 4*tf^3, 3*tf^2, 2*tf, 1, 0;
%               42*tf^5, 30*tf^4, 20*tf^3, 12*tf^2, 6*tf, 2, 0, 0;
%               210*tf^4, 120*tf^3, 60*tf^2, 24*tf, 6, 0, 0, 0];
snap_matrix = [1, t, t^2, t^3, t^4, t^5;
    0, 1, 2*t, 3*t^2, 4*t^3, 5*t^4;
    0, 0, 2, 6*t, 12*t^2, 20*t^3;
    1, tf, tf^2, tf^3, tf^4, tf^5;
    0, 1, 2*tf, 3*tf^2, 4*tf^3, 5*tf^4;
    0, 0, 2, 6*tf, 12*tf^2, 20*tf^3];
state_matrix = [start; zeros(2,3); goal; zeros(2,3)];

% ********************* Compute coefficient matrix ************************
snap_coeffs = snap_matrix \ state_matrix;       
end


%% ********************* Generate pos, vel and acc ************************ 
function [pos, vel, acc] = computeTrajectory(time_matrix, t, coeff_matrix)
% COMPUTE TRAJECTORY
% DO: Compute the desire states of drone according to the time now, and the
%     coefficient matrix above
%     Input:
%       time_matrix: [n, 3] matrix restores the time of each waypoint
%       t: the time now
%       coeff_matrix: [8*(n-1), 3] matrix, the coefficient matrix computed 
%                      above
%
    
% ************************ Compute real state *****************************
[row_num, ~] = find(time_matrix >= t);
coeff = coeff_matrix(6*(row_num(1))-5: 6*(row_num(1)), :);

T_matrix = [1, t, t^2, t^3, t^4, t^5;
             0, 1, 2*t, 3*t^2, 4*t^3, 5*t^4;
              0, 0, 2, 6*t, 12*t^2, 20*t^3];

result_matrix = T_matrix * coeff;
pos = result_matrix(1,:); vel = result_matrix(2,:); acc = result_matrix(3,:);
end
function [pos, vel, acc] = traject_interpolation(t_start, t_finish, t, point_start, point_final)
    pos = zeros(1,3); vel = zeros(1,3); acc = zeros(1,3);
    matrix_Quintic = [1, t_start, t_start^2, t_start^3, t_start^4, t_start^5;
                      0, 1, 2*t_start, 3*t_start^2, 4*t_start^3, 5*t_start^4;
                      0, 0, 2, 6*t_start, 12*t_start^2, 20*t_start^3;
                      1, t_finish, t_finish^2, t_finish^3, t_finish^4, t_finish^5;
                      0, 1, 2*t_finish, 3*t_finish^2, 4*t_finish^3, 5*t_finish^4;
                      0, 0, 2, 6*t_finish, 12*t_finish^2, 20*t_finish^3];
    matrix_Time = [1, t, t^2, t^3, t^4, t^5;
                   0, 1, 2*t, 3*t^2, 4*t^3, 5*t^4;
                   0, 0, 2, 6*t, 12*t^2, 20*t^3;];
    matrix_State = [point_start; zeros(2,3); point_final; zeros(2,3)];
    matrix_Coefficient = matrix_Quintic \ matrix_State;
    matrix_Result = matrix_Time * matrix_Coefficient;
    pos = matrix_Result(1,:)'; vel = matrix_Result(2,:)'; acc = matrix_Result(3,:)';  %pos, vel, acc, 3x1
end

%======NEW-Function to refine the path points with same direction======
% This function is used to refine and delete points with same direction
% within error
function path_new = path_refine_samedirection_new(path, error)  %path nx3 matrix; path_new mx3 matrix
    orie_error = error;
    path_num = size(path, 1);
    path_diff = diff(path);  %calculate the difference between 2 points, (n-1)x3 matrix
    %calculate distance between 2 points
    path_dis = sqrt(path_diff(:, 1).^2 + path_diff(:, 2).^2 + path_diff(:, 3).^2)';
    path_dis = path_dis';  %(n-1)x3 matrix
    path_orie = path_diff ./ path_dis;
    
    path_new = path;
    delete = [];
    for i = 1: size(path_orie)-1
        if path_orie(i,:) > path_orie(i+1,:) - orie_error & path_orie(i,:) < path_orie(i+1,:) + orie_error
            delete = [delete; i+1];
        end
    end
    path_new(delete,:) = [];
end

%======Function to refine the path points which are too close======
% This function is used to delete those points which are too close so that
% to refine path
function path_new = path_refine_closepts(path, dis_min, dis_max)  %path nx3 matrix, dis_min to indicate min distance to filter
    path_new = path;
    path_diff = diff(path); %path_diff (n-1)x3 matrix
    %calculate distance between 2 points
    path_dis = sqrt(path_diff(:, 1).^2 + path_diff(:, 2).^2 + path_diff(:, 3).^2)';
    path_dis = path_dis';  %(n-1)x3 matrix
    path_delete_row = find(path_dis < dis_max & path_dis > dis_min);
    points_delete_row = path_delete_row + 1;
    path_new(points_delete_row, :) = [];
end


