function v = get_particle(data, IDNum)
% 3DFM function  
% Tracking 
% last modified 05/22/05
%  
% Extracts a bead's video tracking data from load_video_tracking.
%  
%  [v] = get_bead(data, IDNum);  
%   
%  where "data" is the output matrix from load_video_tracking
%        "IDNum" is the bead's ID Number 
%   
%  05/09/05 - created; jcribb. How can this guy by found and strangled?
%  05/22/05 - modified to accomodate table or stucture of vectors format
% determine whether the input data is in the table or structure of vectors
% format...
fprintf('entered get_particle\n');
if isfield(data, 'id')
    idx = find(data.id == IDNum);
    idx
    v.id= data.id(idx);
    v.t = data.t(idx);
    v.frame = data.frame(idx);
    v.x = data.x(idx);
    v.y = data.y(idx);
    v.r = data.r(idx);
    if isfield(data,'roll');    v.roll = data.roll(idx);    end;
    if isfield(data,'pitch');   v.pitch= data.pitch(idx);   end;
    if isfield(data,'yaw');     v.yaw  = data.yaw(idx);     end;                    
else
    this_particle = find(data(:,3) == IDNum);
    v = data(this_particle,:);
end
    v
fprintf('leaving get_particle\n)');
return