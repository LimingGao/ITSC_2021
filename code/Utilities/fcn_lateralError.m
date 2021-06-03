%%%%%%%%%%%%%%%%%%%%%  Function fcn_lateralError %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%      find the projection distance on desired trajectory from actual
%      trajectory
%
% Input Variables:
%      Desired_traj    = reference path, the data point should go from start to the end in sequence, format:N*2 matrix (X,Y)
%      Actual_traj     = actual trajectory data, the data point should go from start to the end in sequence,format:N*2 matrix (X,Y)
%                    
% Returned Results:
%      lateral_error  = the lateral offset between two trajectories, format:N*1 array
%
% Example:
%
%
% Processing Flow:
%
% Restrictions/Notes: 
% 
% The following functions are called:
%      none
%
% Author:             Liming Gao
% Created Date:       2020-09-19
% Revisions:
%  
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lateral_error = fcn_lateralError(Desired_traj,Actual_traj)
Md_Desired_trajectory= KDTreeSearcher(Desired_traj); % prepare data tree for knn search
Id_nearest = knnsearch(Md_Desired_trajectory,Actual_traj,'K',2);

Id_nearest = sort(Id_nearest,2); % sort in Ascending Order

% find the projection point on desired path
nearest_point1 = Desired_traj(Id_nearest(:,1),:);
nearest_point2 = Desired_traj(Id_nearest(:,2),:);

% The vector from nearest point1 to nearest point2 on the traversal path
np1_to_np2 = nearest_point2 - nearest_point1;

% The vector from nearest point1 to query point
np1_to_qp = Actual_traj - nearest_point1;

% the length of  vector
np1_to_np2_Length = sum(np1_to_np2.^2,2).^0.5;
np1_to_qp_length = sum(np1_to_qp.^2,2).^0.5;

%pp means proejction point,and np1_to_pp_length is the projection
%length of np1_to_qp to np1_to_np2, i.e. the distance from nearest
%point1 to projection point (the value can be negtive)
np1_to_pp_length = dot(np1_to_qp,np1_to_np2,2)./ np1_to_np2_Length;

np_to_pp_length = sqrt(np1_to_qp_length.^2 - np1_to_pp_length.^2);

lateral_error = np_to_pp_length;

end

