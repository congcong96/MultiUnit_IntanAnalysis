function files = gfn(id, fullpathopt)

% Get file names.

% Convenience function to quickly get a cell array of file names based on
% wildcard calls. Assumes that you're in the directory with the files,
% otherwise calling the full path is necessary. 
% 
%   id: filename patterns, accepts wildcard calls.
% 
% Updated 4/6/17 to include fullpathopt by JS.

if nargin == 1
    fullpathopt = 0;
end

d = dir(id);
files = {d.name}';

if fullpathopt == 1
    files = cellfun(@(x,y) fullfile(y, x), files, {d.folder}', 'UniformOutput', 0);
end

% if length(files) == 1
%     files = files{1};
% end

