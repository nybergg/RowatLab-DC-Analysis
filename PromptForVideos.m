%% PromptForVideos.m
% function [pathNames, videoNames] = PromptForVideos(videoSearchPath)
% Opens a GUI to select videos, user can select a single file at a time and
% ends selection by clicking 'cancel'. Video and path names returned
% are loaded into and returned by the cell arrays 'pathNames' and 'videoNames'.

% Code from Dr. Amy Rowat's Lab, UCLA Department of Integrative Biology and
% Physiology
% Code by Ajay Gopinath (July 2013)

% Inputs
%   - videoSearchPath: specifies the default directory to select videos

% Outputs
%   - pathNames: a string or array of strings that specifies the filepaths
%       of the videos.
%   - videoNames: a string or array of strings that specifies the names of
%       the videos

function [pathNames, videoNames] = PromptForVideos(videoSearchPath)

DEBUG_FLAG = false;

filePath = videoSearchPath;
lastFilename = 1;
pathNames = cell(1,1);
videoNames = cell(1,1);

if(~DEBUG_FLAG)
    idx = 1;
    while (lastFilename ~= 0)
        [fileNames, filePath] = uigetfile('.avi', 'Please select the next video file, if done select cancel.', filePath, 'multiselect', 'on');
        lastFilename = 0;
		if(iscell(fileNames))
			for i = 1:size(fileNames, 2)
				pathNames{idx} = filePath;
				videoNames{idx} = fileNames{i};
				idx = idx + 1;
				lastFilename = fileNames{i};
			end
        elseif(ischar(fileNames))
			pathNames{idx} = filePath;
			videoNames{idx} = fileNames; % fileNames is just a single filename string in this case
			idx = idx + 1;
			lastFilename = fileNames;
		end
    end
else
    % if you want to process all videos in a set of folders, uncomment
    % 'searchPaths' to specify the folders of the videos to process:
    % searchPaths = {'Y:\Kendra\Microfluidics\Raw Data\130618\Mock\7x10um\2psi\';
    %                'Y:\Kendra\Microfluidics\Raw Data\130628\Mock\7x10\2psi\'};
    % if not, manually specify the elements of 'pathNames' and 'videoNames'
    if ~exist('searchPaths', 'var')
        % make sure to delimit path and video names below by semicolons, NOT commas
        pathNames = {'G:\CellVideos'};
        videoNames = {'dev9x10_20X_1200fps_0.6ms_2psi_p9_324_1.avi'};
        
        % go through the specified videos to check if they exist
        for i = 1:length(videoNames)
            if ~(exist(fullfile(pathNames{i}, videoNames{i}), 'file') == 2)
                disp(['Error: ', fullfile(pathNames{i}, videoNames{i}), ' doesnt exist']);
                return;
            end
        end
    else
        idx = 0;
        % go through every path to search
        for i = 1:length(searchPaths)
            pathToSearch = searchPaths{i};
            % throw an error if pathToSearch doesn't exist
            if ~(exist(fullfile(pathToSearch), 'file') == 7)
                disp(['Error: ', pathToSearch, ' doesnt exist']);
                return;
            end
            
            % allVids is the list of all .avi files in pathToSearch
            allVids = dir(fullfile(pathToSearch, '*.avi'));
            
            % store the file/path name of each file in allVids in 'pathNames' 
            % and 'videoNames'
            for j = 1:length(allVids)
                pathNames{j+idx} = pathToSearch;
                videoNames{j+idx} = allVids(j).name;
            end
            idx = idx + length(allVids);
        end
    end
end
