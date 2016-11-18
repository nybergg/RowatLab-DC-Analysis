%% MainCode.m
% Main loop of the cell deformer tracking code.  Designed to input .avi
% video files and track cells through the cell deformer device.  Designed
% to run multiple videos with minimal user interaction.

% Code from Dr. Amy Rowat's Lab, UCLA Department of Integrative Biology and
% Physiology
% Code originally by Bino Varghese (October 2011)
% Updated by David Hoelzle (January 2013)
% Updated by Sam Bruce, Ajay Gopinath, and Mike Scott (July 2013)
% Updated by Kendra Nyberg (May 2014)
% Updated by Kendra Nyberg (April 2015)
% Updated by Kendra Nyberg (June 2015)
% Updated by Kendra Nyberg (November 2016)

% Inputs
%   - .avi files are selected using a GUI
%   - The video names should include, anywhere in the name, the following:
%   1) "devNxM" where N is constriction width / template size to use
%       ex. "dev5x10..." and M is the constriction height       
%   2) "Xfps", where X is the frame rate in frames per second
%       ex. "1200fps..."
%       Example of a properly formatted video name:
%       'dev5x10_1200fps_48hrppb_glass_4psi_20x_0.6ms_p12_041'

% Outputs
%   (Updated: April 2015 by Kendra N)
%   - An excel file with 6 sheets at the specified compiledDataPath
%       1) Total transit time (ms) and unconstricted area (pixels)
%       2) Transit time data (ms)
%       3) Area information at each constriction (pixels)
%       4) Approximate diameter at each constriction (pixels), calculated
%       as the average of major and minor axis of the cell
%       5) Eccentricity of each cell at each constriction
%       6) Retention of each video containing Retention %, # and total #

% Functions called
%   - PromptForVideos   (opens a GUI to load videos)
%   - MakeWaypoints     (Determines the constriction regions)
%   - CellDetection     (Filters the video frames to isolate the cells)
%   - CellTracking      (Labels and tracks cells through the device)
%   - progressbar       (Gives an indication of how much time is left)

% Updated 7/2013 by Mike
%       - Cut out the preprocessing 50 frames (required editing indicies of
%       the call for CellDetection
%       - Rearranged and commented the code to make it clearer
%       - Added the template.  Now MakeWaypoints is automatic and no
%       longer requires defining the cropping and constriction regions
%       - Eliminated redundant inputs and outputs from functions
%       - Eliminated 'segments', nobody used them
% Updated 7/16/2013 by Ajay
%       - separated all logic for prompting/selection of video files to
%       process into function PromptForVideos
%       - improved extraction of frame rates and template sizes from video
%       names by using regular expressions instead of ad-hoc parsing
%       - cleaned up any remaining legacy code and comments
%       - added better output of debugging information
% Updated 3/27/2014 by Kendra
%       - Added retention data calculation
%       - Eliminated discard of transit time data for cells that partially
%       pass through the device. This enables retention calculation and
%       analysis
% Updated 5/4/2014 by Kendra
%       - Added partition of videos by three to increase sample size while
%       keeping variability low.
% Updated 5/6/2014 by Kendra
%       - Added a variable to load filter settings into for agarose
%       detection.
% Updated 4/4/2015 by Kendra
%       - Rearranged excel output to accomadate OriginLab
%       - Added an output with processed raw data (effective cell diameter
%       + c1 transit)
%       - Added lane occlusion matrix - defined as the number of particles
%       in the first constriction of neighboring lanes as a particle 
%       transits through the first constriction.
%       - Additional new output containing:
%             - Effective Diameter
%             - C1 transit
%             - Entry Frame
%             - Exit Frame
%             - Total Occlusion at entry
%             - Total Occlusion at exit
%             - Average Occlusion during transit
%             - Macimum Occlusion during transit


close all
clear variables
clc

%% Initializations
% Allocates an array for the data
numDataCols = 9;
writeData = [];
pairedWriteData = [];
folderNames = [];
compData = [];
pairedData = [];
compiledDataPath = 'C:\Users\Mike\Desktop\';


%% User input dialog box initiation
prompt = {sprintf('1. What type of particles are in your videos?\n(1=Cells, 2=Agrose Gels, or 3=Oil Droplets)'),...
    '2. Do you want to check for correct mask fitting? (Y/N)', ...
    '3. Do you want to use user-defined detection filters? (Y/N)'};
dlg_title = 'User Preferences';
num_lines = 1;
def = {'1','N','N'};
userPref = inputdlg(prompt,dlg_title,num_lines,def); %newid is a function similar to inputdlg and requires download to the toolbox folder


%% Loads video files and prepare any metadata
[pathNames, videoNames] = PromptForVideos('G:\CellVideos\');

% Checks to make sure at least one video was selected for processing
if(isempty(videoNames{1}))
    disp('No videos selected.');
    close all;
    return;
end

%% Extracts the template size and frame rates from the video name.
%   The video names should include, anywhere in the name, the following:
%   1) "devNx" where N is constriction width / template size to use
%       ex. "dev5x10"
%   2) "Mfps", where M is the frame rate in frames per second
%       ex. "1200fps"
% Example of properly formatted video names:
% 'dev5x10_1200fps_48hrppb_glass_4psi_20x_0.6ms_p12_041'

% Allocating matrices
templateSizes = ones(1, length(videoNames));
frameRates = ones(1, length(videoNames));

for i = 1:length(videoNames)
    videoName = videoNames{i};
    if regexp(videoName, 'dev\d*x')
        [j,k] = regexp(videoName, 'dev\d*x');
    elseif regexp(videoName, 'Dev\d*x')
        [j,k] = regexp(videoName, 'Dev\d*x'); % store start/end indices of template size
    else
        disp('File name error.');
        close all;
        return;
    end
    
    [m, n] = regexp(videoName, '\d*fps'); % store start/end indices of frame rate
    templateSize = videoName((j+3):(k-1)); % removes 'dev' at the start, and 'x' at the end
    frameRate = videoName(m:(n-3)); % removes 'fps'  at the end
    
    templateSizes(i) = str2double(templateSize);
    frameRates(i) = str2double(frameRate);
end


%% Loads filter preferences
% filters is a matrix of filter preferences per video each row corresponds
% to the following parameters:
% (1)sample bckgrnd window; (2)backgroundImg averaging pixel range; (3) threshold division--> the greater, the lower the threshold;
% (4)bwareaopen size range cutoff; (5) median filtering radius to remove excess
% noise (originally 5)

if strcmp(userPref(3),'Y') == 1
    disp('Searching for filters.xlsx in first data folder...');
    %Loads excel file for gel detection filter settings
    filters = xlsread([pathNames{1}, 'filters.xlsx']);
    
    % Checks to make sure the filters were located
    if(isempty(filters))
        disp('No filters found.');
        return;
    else
        disp('Filters found.');
    end
    
    % Checks to make sure the filters are the same length as the Video Names
    if(size(filters,2)~=size(videoNames,2))
        disp('The number of filters does not correspond with the number of videos.');
        return;
    end
    
    filters = filters(1:5,:);

elseif strcmp(userPref(1),'2') == 1
    filters = [500; 3; 30; 35; 5];
    filters = repmat(filters, 1, size(videoNames,2));
end


%% Checks for mask alignment
if strcmp(userPref(2),'Y') == 1
    disp('Warning: Videos will not be processed.')
    fprintf('\nAligning masks...\n')
    
    alignment = {'File Path', 'Video Name', 'Y Offset', 'Start Y Coordinate', 'End Y Coordinate'};
    
    for i = 1:length(videoNames)
        % Initializations
        currPathName = pathNames{i};
        currVideoName = videoNames{i};
        currVideo = VideoReader(fullfile(currPathName, currVideoName));
        
        % Calls the MakeWaypoints function to define the constriction region.
        % This function draws a template with a line across each constriction;
        % these lines are used in calculating the transit time
        [mask, lineTemplate, xOffset, yOffset, maskCheck] = MakeWaypoints(currVideo, templateSizes(i), userPref(2));
        
        if i==1
            maskChecks(:,:,:,i) = maskCheck;
        elseif isequal(size(maskCheck),size(maskChecks(:,:,:,1)))
            maskChecks(:,:,:,i) = maskCheck;
        else
            maskChecks(:,:,:,i) = zeros(size(maskChecks(:,:,:,1)));
        end
        
        startPos = 33 + yOffset;
        endPos = 65 + yOffset;
        
        alignment = vertcat(alignment, {currPathName, currVideoName, yOffset, startPos, endPos});
        
    end
    
    % xlswrite(fullfile(pathNames{1}, 'Manual Transit Time Positions'), alignment, 'Sheet1', 'A1');
    
    maskChecksMov = immovie(maskChecks);
    implay(maskChecksMov,4)
    set(findall(0, 'tag', 'spcui_scope_framework'), 'position', [500 300 600 400]);
    fprintf('Done.\n')
    return
end


%% Create the folder in which to store the output data
% The output folder name is a subfolder in the folder where the first videos
% were selected. The folder name contains the time at which processing is
% started.

outputFolderName = fullfile(pathNames{1}, ['processed_', datestr(now, 'mm-dd-YY_HH-MM')]);
if ~(exist(outputFolderName, 'file') == 7)
    mkdir(outputFolderName);
end

lastPathName = pathNames{i};


%% Initializes a progress bar
progressbar('Overall', 'Cell detection', 'Cell tracking');
tStart = tic;


%% Iterates through videos to filter, analyze, and output the compiled data
for i = 1:length(videoNames)
    % Initializations
    currPathName = pathNames{i};
    outputFilename = fullfile(outputFolderName, regexprep(currPathName, '[^a-zA-Z_0-9-]', '~', 'all'));
    currVideoName = videoNames{i};
    currVideo = VideoReader(fullfile(currPathName, currVideoName));
    disp(['==Video ', num2str(i), '==']);
 
    startFrame = 1;
    endFrame = currVideo.NumberOfFrames;
    
    % Calls the MakeWaypoints function to define the constriction region.
    % This function draws a template with a line across each constriction;
    % these lines are used in calculating the transit time
    [mask, lineTemplate, xOffset] = MakeWaypoints(currVideo, templateSizes(i));

    
    % Calls CellDetection to filter the images and store them in
    % 'processedFrames'.  These stored image are binary and should
    % (hopefully) only have the cells in them
    if strcmp(userPref(1),'1') == 1 % Videos containing cells
        [processedFrames] = CellDetection(currVideo, startFrame, endFrame, currVideoName, mask);
    elseif strcmp(userPref(1),'2') == 1 % Videos containing gels
        [processedFrames] = GelDetection(currVideo, startFrame, endFrame, currVideoName, mask, filters);
    elseif strcmp(userPref(1),'3') == 1 % Videos containing oil droplets
        [processedFrames] = OilDetection(currVideo, startFrame, endFrame, currVideoName, mask);
    end  
        
    % Calls CellTracking to track the detected cells.
    [occl, comp, debuglaneData, cellInfo, paired] = CellTracking(startFrame, endFrame, frameRates(i), lineTemplate, processedFrames, xOffset, currVideoName);
    progressbar((i/(size(videoNames,2))), 0, 0)
    
    
    % If data is generated (cells are found and tracked through the device)
    
    if (~isempty(comp))
        compData = vertcat(compData, comp);
    end
    
    if (~isempty(paired))
        pairedData = vertcat(pairedData, paired);
    end
    
    
    lastPathName = currPathName;
    
    %% Prepare for data output
    if((i == length(videoNames)) || ~strcmp(pathNames{i}, pathNames{i+1}))
        
        % Prepare data for Compiled Raw Data
        [writeData, folderNames] = writeCOMPoutput(pathNames, writeData, i,folderNames, compData);
        compData = [];
        
        [pairedWriteData, folderNames] = writeCOMPoutput(pathNames, pairedWriteData, i,folderNames, pairedData);
        pairedData = [];
    
    end
end

%Output excel file ready for Compiled Raw Data
colHeader3 = {'Diameter (um)', 'C1 Transit Time (ms)', 'Entry Frame', 'Exit Frame', 'Entry Occlusion', 'Exit Occlusion', 'Average Occlusion', 'Max Occlusion', 'Lane', 'VideoName','filepath:'};
compOutput = fullfile(outputFolderName, 'Compiled Processed Data');
xlswrite(compOutput,colHeader3,'Sheet1','A1');
xlswrite(compOutput,writeData,'Sheet1','A2');


%% Output debugging information
totalTime = toc(tStart);
avgTimePerVideo = totalTime/length(videoNames);

fprintf('\n\n===========\n');
disp(['Total time to analyze ', num2str(length(videoNames)), ' video(s): ', num2str(totalTime), ' secs']);
disp(['Average time per video: ', num2str(avgTimePerVideo), ' secs']);
fprintf('\n\nOutputting metadata...\n');

runOutputPaths = unique(pathNames);
for i = 1:length(runOutputPaths)
    runOutputFile = fopen(fullfile(runOutputPaths{i}, 'process_log.txt'), 'wt+', 'l', 'UTF-8');
    vidIndices = strcmp(runOutputPaths{i}, pathNames);
    vidsProcessed = videoNames(vidIndices);
    
    fprintf(runOutputFile, '%s\n\n', 'The following files were processed from this folder:');
    fprintf(runOutputFile, '%s\n', '============');
    for j = 1:length(vidsProcessed)
        fprintf(runOutputFile, '%s\n', vidsProcessed{j});
    end
    fprintf(runOutputFile, '%s\n\n', '============');
    
    fprintf(runOutputFile, '%s%s\n', 'Processing was finished at: ', datestr(now, 'mm-dd-YY HH:MM:SS'));
    fprintf(runOutputFile, '%s%s\n', 'Output files are located at: ', outputFolderName);
    
    fclose(runOutputFile);
end


% %add occl probability
% occl(1:8470,8) = 1;
% occl(1:10141,12) = 1;
% occl(:,17) = sum(occl(:,1:16),2);
occl(:,18) = (16-occl(:,17))/16;
occl(:,19) = (occl(:,17))/16;
occl(:,20) = 5*(1:length(occl));

disp('Done.')
