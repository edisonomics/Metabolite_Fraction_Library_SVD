function setup()
% SETUP Add all repository folders (excluding hidden directories) to the MATLAB path.
%   Run once from the repository root before executing any workflow script.

if ~isfolder('src')
    error('setup:WrongDirectory', ...
        'Run setup from the repository root directory (the folder containing src/).');
end

repo_root = pwd;
all_paths = strsplit(genpath(repo_root), pathsep);
all_paths = all_paths(~cellfun('isempty', all_paths));

% Exclude hidden directories (.git, .github, .vscode, etc.) and archive/
rel_paths = erase(all_paths, repo_root);
keep      = ~contains(rel_paths, [filesep '.']) & ...
            ~startsWith(rel_paths, [filesep 'archive']);

addpath(strjoin(all_paths(keep), pathsep));
fprintf('Repository added to MATLAB path: %s\n', repo_root);
end
