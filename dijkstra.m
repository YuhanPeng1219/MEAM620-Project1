function [path, num_expanded] = dijkstra(map, start, goal, astar)
% DIJKSTRA Find the shortest path from start to goal.
%   PATH = DIJKSTRA(map, start, goal) returns an mx3 matrix, where each row
%   consists of the (x, y, z) coordinates of a point on the path. The first
%   row is start and the last row is goal. If no path is found, PATH is a
%   0x3 matrix. Consecutive points in PATH should not be farther apart than
%   neighboring voxels in the map (e.g. if 5 consecutive points in PATH are
%   co-linear, don't simplify PATH by removing the 3 intermediate points).
%
%   PATH = DIJKSTRA(map, start, goal, astar) finds the path using euclidean
%   distance to goal as a heuristic if astar is true.
%
%   [PATH, NUM_EXPANDED] = DIJKSTRA(...) returns the path as well as
%   the number of nodes that were expanded while performing the search.
%   
% paramaters:
%   map     - the map object to plan in
%   start   - 1x3 vector of the starting coordinates [x,y,z]
%   goal:   - 1x3 vector of the goal coordinates [x,y,z]
%   astar   - boolean use astar or dijkstra

tic;
%% initial program 
% if Astar is default, set false
if nargin < 4
    astar = false;
end

% initial output
path = [];
num_expanded = 0;

%% Check whether start and goal points are in boundary
if boundary_check([start; goal], map)
    warning('The initial points are out of boundary!!!');
    return
end

%% Locate start and goal points
start_ind = pos2ind(map, start);
goal_ind = pos2ind(map, goal);
goal_sub = pos2sub(map, goal);
start_sub = pos2sub(map, start);

%% Check whether start and goal points are valid
% The start and gaol points can not collide with obstacles
input_ind = [start_ind; goal_ind];
if any(map.occgrid(input_ind) == 1) 
   warning('Start or goal point collide with blocks, check first!!') ;
   return 
end


%% Prepare for path finding
% Get the information of cells in the map grid.
cell_total = numel(map.occgrid);
[cell_row, cell_col, cell_hei] = size(map.occgrid);

% Set all the cost of cell infinity
node_cost = inf(cell_row, cell_col, cell_hei);
Astar_nodecost = inf(cell_row, cell_col, cell_hei);

% initial the parents cell and current cell
Que = ones(cell_total,1);
parents = zeros(cell_total, 1);
expand_node = zeros(cell_total, 1);
current_ind = start_ind;
node_cost(current_ind) = 0;
Astar_nodecost(current_ind) = Hueristic(start_sub, goal_sub, map);

%% Main loop to extend the search tree
if astar == false
    while true
        % Find the cell in the cost graphy with minimum distance to start point
        % Mark the index 0 in the Que matrix
        % Determine whether the goal is found
        [min_cost, current_ind] = min(node_cost(:));
        if current_ind == goal_ind
            break
        end
        Que(current_ind) = 0;
        expand_node(current_ind) = 1;
        
        % If no path is found
        if isinf(min_cost)
            warning('No path exists!!')
            return
        end
        
        % convert the index of current node into subscripts
        [i,j,k] = ind2sub(size(node_cost), current_ind);
        current_sub = [i, j, k];
        node_cost(current_ind) = inf;
        
        % find the neighbors
        children_sub = find_children(i, j, k);
        
        % check whether the subscripts are in map bound
        isbelow = children_sub < 1;
        isover = children_sub > [cell_row, cell_col, cell_hei];
        outrange = isbelow + isover;
        
        % omit those are out of range
        [omit_row, ~] = find(outrange > 0);
        children_sub(omit_row,:) = [];
        
        % Get index of children
        children_ind = sub2ind(size(map.occgrid),children_sub(:,1),...
                        children_sub(:,2), children_sub(:,3));

         % check whether the subscripts are in obstacle
        Check_result = map.occgrid (children_ind);
        [omit_collision, ~] = find(Check_result > 0);
        
        % check whether has expanded
        expanded = expand_node(children_ind);
        [omit_expand, ~] = find(expanded > 0);
        
        % cancel rows in obstacle
        omit = [omit_expand; omit_collision];
        children_sub(omit,:) = [];
        
        children_ind = sub2ind(size(map.occgrid),children_sub(:,1),...
                        children_sub(:,2), children_sub(:,3));
                    
        % Do: collision detect(using cell subscripts)
        % Debug: plot children
        % Compute the cost for children cell
        % Mark new parents
        % Count for expended nodes                    
        new_cost = min_cost + distToParents(children_sub, current_sub, map);
        [should_change,~] = find(node_cost(children_ind) > new_cost);
        node_cost(children_ind(should_change)) = new_cost(should_change);
        parents(children_ind(should_change)) = current_ind;
        expand_node(children_ind(should_change)) = 1;

%         for n = 1 : size(children_ind)
%             if collisionCell_check (children_sub(n, :), map) ~= 1 && ...
%                     children_ind(n) ~= start_ind && Que(children_ind(n)) ~= 0
%                 % DEBUG here !!
% %                          hold on; plot3(children_sub(n, 1), children_sub(n, 2), children_sub(n, 3), 'bo');
% %                          drawnow;
%                 new_cost = min_cost + distToParents(children_sub(n, :), current_sub, map);
%                 if node_cost(children_ind(n)) > new_cost
%                     node_cost(children_ind(n)) = new_cost;
%                     parents(children_ind(n)) = current_ind;
%                     expand_node(children_ind(n)) = 1;
%                 end
%             end
%         end
    end
else
    while true
        % Find the cell in the cost graphy with minimum distance to start point
        % Mark the index 0 in the Que matrix
        % Determine whether the goal is found
        [~, current_ind] = min(Astar_nodecost(:));
        min_cost = node_cost(current_ind);
        if current_ind == goal_ind
            break
        end
        Que(current_ind) = 0;
        expand_node(current_ind) = 1;
        
        % If no path is found
        if isinf(min_cost)
            warning('No path exists!!')
            return
        end
        
        % convert the index of current node into subscripts
        % Change current node to inf
        [i,j,k] = ind2sub(size(Astar_nodecost), current_ind);
        current_sub = [i, j, k];
        Astar_nodecost(current_ind) = inf;
        
        % find the neighbors
        children_sub = find_children(i, j, k);
       
        % check whether the subscripts are in map bound
        isbelow = children_sub < 1;
        isover = children_sub > [cell_row, cell_col, cell_hei];
        outrange = isbelow + isover;
       
        % omit those are out of range
        [omit_row, ~] = find(outrange > 0);
        children_sub(omit_row,:) = [];
        
        % Get index of children
        children_ind = sub2ind(size(map.occgrid),children_sub(:,1),...
                        children_sub(:,2), children_sub(:,3));
        
        % check whether the subscripts are in obstacle
        Check_result = map.occgrid (children_ind);
        [omit_collision, ~] = find(Check_result > 0);
        
        % check whether has expanded
        expanded = expand_node(children_ind);
        [omit_expand, ~] = find(expanded > 0);
        
        % cancel rows in obstacle
        omit = [omit_expand; omit_collision];
        children_sub(omit,:) = [];
        
        children_ind = sub2ind(size(map.occgrid),children_sub(:,1),...
                        children_sub(:,2), children_sub(:,3));
                    
        % Do: collision detect(using cell subscripts)
        % Debug: plot children
        % Compute the cost for children cell
        % Mark new parents
        % Count for expended nodes
%         for n = 1 : size(children_ind)
%             if map.occgrid(children_ind) ~= 1 && ...
%                     children_ind(n) ~= start_ind && expand_node(children_ind(n)) ~= 1
%                 % DEBUG here !!
% %                    hold on; 
% %                    pos_children = sub2pos(map, children_sub(n, :));
% %                    plot3(pos_children( 1), pos_children(2),pos_children(3), 'bo');
% %                    drawnow;
%                 new_cost = min_cost + distToParents(children_sub(n, :), current_sub, map);
%                 Astar_new = new_cost + Hueristic(children_sub(n, :), goal_sub, map);
%                 if node_cost(children_ind(n)) > new_cost
%                     node_cost(children_ind(n)) =  new_cost;
%                     Astar_nodecost(children_ind(n)) = Astar_new;
%                     parents(children_ind(n)) = current_ind;
%                     expand_node(children_ind(n)) = 1;
%                 end
%             end
%         end

        new_cost = min_cost + distToParents(children_sub, current_sub, map);
        Astar_new = new_cost + Hueristic(children_sub, goal_sub, map);
        [should_change,~] = find(Astar_nodecost(children_ind) > Astar_new);
        node_cost(children_ind(should_change)) = new_cost(should_change);
        Astar_nodecost(children_ind(should_change)) = Astar_new(should_change);
        parents(children_ind(should_change)) = current_ind;
        expand_node(children_ind(should_change)) = 1;
%         if node_cost(children_ind(n)) > new_cost
%             node_cost(children_ind(n)) =  new_cost;
%             Astar_nodecost(children_ind(n)) = Astar_new;
%             parents(children_ind(n)) = current_ind;
%             expand_node(children_ind(n)) = 1;
%         end
    end
end
%% Generate path according to parents
path_ind = [];
if (isinf(node_cost(goal_ind)))
    warning('No path is found!!');
    return
else
    path_ind = [path_ind; goal_ind];
    while (parents(path_ind(1)) ~= 0)
        path_ind = [parents(path_ind(1)); path_ind];
    end
    [pathsub_x, pathsub_y, pathsub_z] = ind2sub(size(map.occgrid), path_ind);
    pathsub = [pathsub_x, pathsub_y, pathsub_z];
    path = sub2pos(map, pathsub);
end

path(1,:) = start;
path(end, :) = goal;
%% Compute total expanded nodes
num_expanded = sum(expand_node);
toc;
end

% %% Function used to check point collision
% function result = collisionPoint_check (points, map)
%    for t = 1 : size(map.blocks)
%       obs = map.blocks(t, 1:6);
%       obs(1:3) = obs(1:3) - map.margin;
%       obs(4:6) = obs(4:6) + map.margin;
%       check_x = points(:, 1) >  obs(1) & points(:,1) < obs(4);
%       check_y = points(:, 2) >  obs(2) & points(:,2) < obs(5);
%       check_z = points(:, 3) >  obs(3) & points(:,3) < obs(6);
%       
%       if any(all([check_x, check_y, check_z],2)) == 1
%           result = 1;
%           break;
%       else
%           result = 0;
%       end
%    end
% end
% 
% %% Function used to check cell collision
% function result = collisionCell_check (cell_sub, map)
% result = zeros(size(cell_sub,1), 1);
% check_x = zeros(size(cell_sub,1), size(map.blocks,1));
% check_y = zeros(size(cell_sub,1), size(map.blocks,1));
% check_y = zeros(size(cell_sub,1), size(map.blocks,1));
%    for t = 1 : size(map.blocks)
%       obs = map.blocks(t, 1:6);
%       obs(1:3) = obs(1:3) - map.margin;
%       obs(4:6) = obs(4:6) + map.margin;
%       pos = sub2pos(map, cell_sub);
%       check_x(:,t) = pos(:, 1) >  obs(1) & pos(:,1) < obs(4);
%       check_y(:,t) = pos(:, 2) >  obs(2) & pos(:,2) < obs(5);
%       check_z(:,t) = pos(:, 3) >  obs(3) & pos(:,3) < obs(6);     
% %       if all([check_x, check_y, check_z])== 1
% %           result = 1;
% %           break;
% %       else
% %           result = 0;
% %       end
%    end
%     result = all(check_x & check_y & check_z);
% %    pos = sub2pos(map, cell_sub);
% end

%% Function to check boundry
 function result = boundary_check(input_points, map)
    if all(all(input_points >= map.bound_xyz(1:3) == 1)) && ...
       all(all(input_points <= map. bound_xyz(4:6) == 1))
       result = 0 ;
    else
       result = 1;
    end
end
%% Function to find children according to parent
 function children = find_children(i, j, k)
    children = [i - 1, j - 1, k - 1;
        i - 1, j - 1,     k;
        i - 1, j - 1, k + 1;
        i - 1,     j, k - 1;
        i - 1,     j,     k;
        i - 1,     j, k + 1;
        i - 1, j + 1, k - 1;
        i - 1, j + 1,     k;
        i - 1, j + 1, k + 1;
        i, j - 1, k - 1;
        i, j - 1,     k;
        i, j - 1, k + 1;
        i,     j, k - 1;
        i,     j, k + 1;
        i, j + 1, k - 1;
        i, j + 1,     k;
        i, j + 1, k + 1;
        i + 1, j - 1, k - 1;
        i + 1, j - 1,     k;
        i + 1, j - 1, k + 1;
        i + 1,     j, k - 1;
        i + 1,     j,     k;
        i + 1,     j, k + 1;
        i + 1, j + 1, k - 1;
        i + 1, j + 1,     k;
        i + 1, j + 1, k + 1];

% children = [ i - 2, j - 2, k - 2;
%  i - 2, j - 2, k - 1;
%  i - 2, j - 2,     k;
%  i - 2, j - 2, k + 1;
%  i - 2, j - 2, k + 2;
%  i - 2, j - 1, k - 2;
%  i - 2, j - 1, k - 1;
%  i - 2, j - 1,     k;
%  i - 2, j - 1, k + 1;
%  i - 2, j - 1, k + 2;
%  i - 2,     j, k - 2;
%  i - 2,     j, k - 1;
%  i - 2,     j,     k;
%  i - 2,     j, k + 1;
%  i - 2,     j, k + 2;
%  i - 2, j + 1, k - 2;
%  i - 2, j + 1, k - 1;
%  i - 2, j + 1,     k;
%  i - 2, j + 1, k + 1;
%  i - 2, j + 1, k + 2;
%  i - 2, j + 2, k - 2;
%  i - 2, j + 2, k - 1;
%  i - 2, j + 2,     k;
%  i - 2, j + 2, k + 1;
%  i - 2, j + 2, k + 2;
%  i - 1, j - 2, k - 2;
%  i - 1, j - 2, k - 1;
%  i - 1, j - 2,     k;
%  i - 1, j - 2, k + 1;
%  i - 1, j - 2, k + 2;
%  i - 1, j - 1, k - 2;
%  i - 1, j - 1, k - 1;
%  i - 1, j - 1,     k;
%  i - 1, j - 1, k + 1;
%  i - 1, j - 1, k + 2;
%  i - 1,     j, k - 2;
%  i - 1,     j, k - 1;
%  i - 1,     j,     k;
%  i - 1,     j, k + 1;
%  i - 1,     j, k + 2;
%  i - 1, j + 1, k - 2;
%  i - 1, j + 1, k - 1;
%  i - 1, j + 1,     k;
%  i - 1, j + 1, k + 1;
%  i - 1, j + 1, k + 2;
%  i - 1, j + 2, k - 2;
%  i - 1, j + 2, k - 1;
%  i - 1, j + 2,     k;
%  i - 1, j + 2, k + 1;
%  i - 1, j + 2, k + 2;
%      i, j - 2, k - 2;
%      i, j - 2, k - 1;
%      i, j - 2,     k;
%      i, j - 2, k + 1;
%      i, j - 2, k + 2;
%      i, j - 1, k - 2;
%      i, j - 1, k - 1;
%      i, j - 1,     k;
%      i, j - 1, k + 1;
%      i, j - 1, k + 2;
%      i,     j, k - 2;
%      i,     j, k - 1;
%      i,     j,     k;
%      i,     j, k + 1;
%      i,     j, k + 2;
%      i, j + 1, k - 2;
%      i, j + 1, k - 1;
%      i, j + 1,     k;
%      i, j + 1, k + 1;
%      i, j + 1, k + 2;
%      i, j + 2, k - 2;
%     i, j + 2, k - 1;
%      i, j + 2,     k;
%      i, j + 2, k + 1;
%      i, j + 2, k + 2;
%  i + 1, j - 2, k - 2;
%  i + 1, j - 2, k - 1;
%  i + 1, j - 2,     k;
%  i + 1, j - 2, k + 1;
%  i + 1, j - 2, k + 2;
%  i + 1, j - 1, k - 2;
%  i + 1, j - 1, k - 1;
%  i + 1, j - 1,     k;
%  i + 1, j - 1, k + 1;
%  i + 1, j - 1, k + 2;
%  i + 1,     j, k - 2;
%  i + 1,     j, k - 1;
%  i + 1,     j,     k;
%  i + 1,     j, k + 1;
%  i + 1,     j, k + 2;
%  i + 1, j + 1, k - 2;
%  i + 1, j + 1, k - 1;
%  i + 1, j + 1,     k;
%  i + 1, j + 1, k + 1;
%  i + 1, j + 1, k + 2;
%  i + 1, j + 2, k - 2;
%  i + 1, j + 2, k - 1;
%  i + 1, j + 2,     k;
%  i + 1, j + 2, k + 1;
%  i + 1, j + 2, k + 2;
%  i + 2, j - 2, k - 2;
%  i + 2, j - 2, k - 1;
%  i + 2, j - 2,     k;
%  i + 2, j - 2, k + 1;
%  i + 2, j - 2, k + 2;
%  i + 2, j - 1, k - 2;
%  i + 2, j - 1, k - 1;
%  i + 2, j - 1,     k;
%  i + 2, j - 1, k + 1;
%  i + 2, j - 1, k + 2;
%  i + 2,     j, k - 2;
%  i + 2,     j, k - 1;
%  i + 2,     j,     k;
%  i + 2,     j, k + 1;
%  i + 2,     j, k + 2;
%  i + 2, j + 1, k - 2;
%  i + 2, j + 1, k - 1;
%  i + 2, j + 1,     k;
%  i + 2, j + 1, k + 1;
%  i + 2, j + 1, k + 2;
%  i + 2, j + 2, k - 2;
%  i + 2, j + 2, k - 1;
%  i + 2, j + 2,     k;
%  i + 2, j + 2, k + 1;
%  i + 2, j + 2, k + 2];

% children = [
%         i - 1, j - 1,     k;
%         i - 1,     j,     k;
%         i - 1, j + 1,     k;
%         i, j - 1,     k;
%         i,     j, k - 1;
%         i,     j, k + 1;
%         i, j + 1,     k;
%         i + 1, j - 1,     k;
%         i + 1,     j,     k;
%         i + 1, j + 1,     k;
%     ];
end

%% Function to compute distance between parent and children
function dist = distToParents(children, parent, map)
   diff_cell = abs(children - parent);
%    dist = sqrt ((map.res_xyz(1) * abs(diff_cell(1)))^2 + (map.res_xyz(2) * ...
%             abs(diff_cell(2)))^2 + (map.res_xyz(3) * abs(diff_cell(3)))^2);
%     dist = map.res_xyz(1) * abs(diff_cell(1))+(map.res_xyz(2) * ...
%              abs(diff_cell(2)))+ (map.res_xyz(3) * abs(diff_cell(3)));
   dist = sqrt ((map.res_xyz(1) * diff_cell(:,1)).^2 + (map.res_xyz(2) * ...
            (diff_cell(:,2))).^2 + (map.res_xyz(3) * diff_cell(:,3)).^2);

end

%% Function to compute Huristics distance
function dist = Hueristic(children, goal, map)
     diff_cell = abs(goal - children);
    diff_pos = [map.res_xyz(1) * diff_cell(:,1), map.res_xyz(2)...
                    * diff_cell(:,2), map.res_xyz(3) * diff_cell(:,3)];
%     [dis, ind] = min(diff_pos,1);
%     diff_pos(:, ind) = [];
%     [dis2, ind2] = min(diff_pos,1);
%     diff_pos(:, ind2) = [];
%     dist = sqrt(3*dis.^2) + sqrt(2*(dis2 - dis).^2) + (diff_pos - dis - dis2);
%         dist = sum(diff_pos, 2);
%            dist = sqrt(diff_pos(:,1).^2+ diff_pos(:,2).^2 + diff_pos(:,3));
          dist = sqrt(diff_pos(:,1).^2+ diff_pos(:,2).^2 + diff_pos(:,3).^2);
end
