%% CellTracking.m
% function [transitData] = CellTracking(numFrames, framerate, template, processedFrames, xOffset)
% Inputs the processed frames of a video, labels them, and stores the data.
% Now evaluates data before storing it, eliminating storage of data that is
% later unused.  Calls ProcessTrackingData to process the raw data.

% Code from Dr. Amy Rowat's Lab, UCLA Department of Integrative Biology and
% Physiology
% Code originally by Bino Varghese (AnalysisCodeBAV.m) (October 2011)
% Updated by David Hoelzle (January 2013)
% Updated by Sam Bruce, Ajay Gopinath, and Mike Scott (July 2013)

% Inputs
%   - numFrames: an integer, the number of frames in the video
%   - framerate: an integer, the framerate of the video
%   - template: a logical array which specifies the horizontal lines used
%       in tracking the cells
%   - processedFrames: a 3 dimensional array that stores the processed
%       frames
%   - xOffset: an integer that stores the number of pixels the template is
%       offset, used to calculate the x-coordinates of the lanes

% Outputs
%   - transitData: an array of data with dimensions (n x 8 x 4) where n
%   is the number of cells found in the video
%       - transitData(:,:,1) is the transit time data
%       - transitData(:,:,2) is the area data
%       - transitData(:,:,3) is the equivalent diameter data
%       - transitData(:,:,4) is the eccentricity data

% Functions called
%   - ProcessTrackingData   (processes the raw data to track the cells)

function [occl, comp, debuglaneData, cellInfo, paired] = CellTracking(startFrame, endFrame, framerate, template, processedFrames, xOffset, currVideoName)
%% Testing
% processedFrames = load('testFrames.mat');
% processedFrames = processedFrames.processedFrames;
% numFrames = (length(processedFrames));
% xOffset = 5;
% framerate = 200;
% template = load('lineTemplate.mat');
% template = template.lineTemplate;

%%

% progressbar([],[],0)
% Change WRITEVIDEO_FLAG to true in order to print a video of the output,
% defaults to false.
WRITEVIDEO_FLAG = false;

% %% Debugging Initializations
% folderName = 'Y:\Kendra\Oil Droplets\131218 - Silicone Oil with Varying Viscosities I\10 cSt\';
% videoName = 'dev5x5_400fps_1hppt_1o10PDMS_4psi_20x_0.4ms__35umfilter_002.avi';
% 
% cellVideo = VideoReader([folderName, videoName]);
% startFrame = 1;
% endFrame = cellVideo.NumberOfFrames;
% 
% % to determine mask variable
% [j,k] = regexp(videoName, 'dev\d*x'); % store start/end indices of template size
% templateSize = videoName((j+3):(k-1)); % removes 'dev' at the start, and 'x' at the end
% userPref = {'N','N','N', '1', 'N'};
% [mask, template, xOffset] = MakeWaypoints(cellVideo, templateSize, userPref(5));
% 
% load('processedFrames.mat', 'processedFrames');
% 
% framerate = str2double('400');



%% Initializations
counter = 0;
line = 0;
write = true;
numFrames = endFrame-startFrame+1;
% HARD CODED x coordinates for the center of each lane (1-16), shifted by
% the offset found in 'MakeWaypoints'
laneCoords = [16 48 81 113 146 178 210 243 276 308 341 373 406 438 471 503] + xOffset;

% The cell structure 'cellInfo' is an important structure that stores
% information about each cell.  It contains 16 arrays, one for each lane in
% the device.  Each array is initialized to a default length
% of 1,500 rows, and the code checks that the array is not full at the end
% of each loop.  If the array is full, it is enlarged. The columns are:
%   1) Frame number
%   2) Cell label number
%   3) Grid line that the cell intersects
%   4) Cell area (in pixels)
%   5) Major axis length
%   6) Minor axis length
cellInfo = cell(1,16);
for ii = 1:16
    cellInfo{ii} = zeros(300,6);
end

% In order to remember which index to write to in each of the arrays in
% cellInfo, a counter variable is needed.  laneIndex gives the index for
% each lane.
laneIndex = ones(1,16);

% The array checkingArray is 'number of horizontal lines' x 'number of
% lanes'.  Each time a cell is found, the position (lane and line) is
% known.  These are used as indicies (for instance, a cell at line 2 in
% lane 4 will check checkingArray(2,4)).  At each position, the last frame
% at which a cell was previously found in that position is stored.  The
% cell is only stored as a cell if no cell was found in sequential previous
% frames.  If the frame stored at that position is 0 or 2 less than the
% current frame, the cell is counted (write is turned true)
checkingArray = zeros(8,16);

%% Labels each grid line in the template from 1-8 starting at the top
[tempmask, ~] = bwlabel(template);

% Preallocates an array to store the y coordinate of each line
lineCoordinate = zeros(1,8);

% Uses the labeled template to find the y coordinate of each line
for jj = 1:8
    q = regionprops(ismember(tempmask, jj), 'PixelList');
    lineCoordinate(jj) = q(1,1).PixelList(1,2);
end
clear tempmask;

%% Opens a videowriter object if needed
if(WRITEVIDEO_FLAG)
    outputVideo = VideoWriter('Y:\Kendra\Agarose Microgels\140226 - Agarose 1.5 1.75 and 2.0 percent at 8psi\1.5 percent\output_video.avi');
    outputVideo.FrameRate = framerate;
    open(outputVideo)
end

%% Cell Labeling
% This loop goes through the video frame by frame and labels all of the
% cells.  It stores (in cellInfo), the centroids and line intersection of
% each cell.
for ii = startFrame:endFrame
    % currentFrame stores the frame that is currently being processed
    currentFrame = processedFrames(:,:,ii);
    % Allocates a working frame (all black).  Any cell in the current
    % frame that is valid (touching a line and of a certain size) will be
    % added into working frame.
    workingFrame = false(size(currentFrame));
    
    % If the current frame has any objects in it.  Skips any empty frames.
    if any(currentFrame(:) ~= 0)
        %% Label the current frame
        % Count number of cells in the frame and label them
        % (numLabels gives the number of cells found in that frame)
        [labeledFrame, numLabels] = bwlabel(currentFrame);
        % Compute their centroids
        cellCentroids = regionprops(labeledFrame, 'centroid', 'area', 'MajorAxisLength', 'MinorAxisLength');
        
        %% Check which line the object intersects with
        for jj = 1:numLabels
            currentRegion = ismember(labeledFrame, jj);
            % Goes through each labeled region in the current frame, and finds
            % the line that the object intersects with
            % Goes through each labeled region in the current frame, and finds
            % the line that the object intersects with
            % The reason for reverse iteration through the lines is
            % because a cell passaging may very well intersect both lines 1
            % (the line for finding "unconstricted cell size") and line 2,
            % the cell may very well intersect both at the same time.
            % To fix this, we want to find when the cell is intersecting
            % line 2 no matter whether it is also intersecting line 1,
            % and iterating in reverse will let us check for line 2's
            % intersection before that of line 1.
            
            linesOneAndTwo = false;
            for line = 8:-1:1
                % Find their intersection
                if(sum(currentRegion(lineCoordinate(line),:)) ~= 0)
                    counter = counter + 1;
                    % Makes sure if the cell is touching line 2, it is not also
                    % touching line 1 (often line 1 is skipped due to its
                    % proximity to line 2)
                    if(line == 2)
                        if(sum(currentRegion(lineCoordinate(1),:)) ~= 0)
                            linesOneAndTwo = true;
                        end
                    end
                    % Breaks to preserve line, the line intersection
                    write = true;
                    break;
                end
                
                
                % If the cell is not touching any of the lines, set
                % line = 0, so it is not included in the array
                % cellInfo
                if(ii == 1)
                    if(line == 1)
                        write = true;
                        break;
                    end
                else
                    if(line == 1)
                    write = false;
                        break;
                    end
                end
            end
            
            if(counter > 0 && write == true && line ~= 0)
                % Determines which lane the current cell is in
                [offCenter, lane] = min(abs(laneCoords-cellCentroids(jj,1).Centroid(1)));
                
                % Now that line and lane are both known, checks the array
                % 'checkingArray' to see if the cell should be stored.
                % There are two possibilities:
                %       1) The element of 'checking array' contains the
                %       previous frame number. In this case, the frame
                %       value in 'checkingArray' is updated, but the cell
                %       is not stored.
                %       2) The element of 'checking array' does not contain
                %       the previous frame number, or contains zero.  In
                %       this case, the frame value is stored in 'checking
                %       array' and the cell is stored in the appropriate
                %       array in cellInfo.
                % In case 1:
                if(offCenter <= 12)
                    if(checkingArray(line,lane) == ii - 1 && (checkingArray(line,lane) ~= 0 || ii ~= 1))
                        % If the cell was in the same place as the line before
                        % And a cell was previously found at this line
                        % ~=0 since if the cell is the first to be found on
                        % line 1 in that lane, it will be zero (if in frame 1)
                        checkingArray(line,lane) = ii;
                    else
                        % If the cell touches lines one and two at the same
                        % time, save information at line 1 and line 2!
                        if(linesOneAndTwo)
                            % Frame number
                            cellInfo{lane}(laneIndex(lane),1) = ii;
                            % Cell number
                            cellInfo{lane}(laneIndex(lane),2) = counter;
                            % Line intersection
                            cellInfo{lane}(laneIndex(lane),3) = 1;
                            % Saves the area of the cell in pixels
                            cellInfo{lane}(laneIndex(lane),4) = cellCentroids(jj,1).Area(1);
                            % Store the length of the major axis
                            cellInfo{lane}(laneIndex(lane),5) = cellCentroids(jj,1).MajorAxisLength;
                            % Store the length of the minor axis
                            cellInfo{lane}(laneIndex(lane),6) = cellCentroids(jj,1).MinorAxisLength;
                            
                            % Updates the checking array and lane index
                            checkingArray(1,lane) = ii;
                            laneIndex(lane) = laneIndex(lane) + 1;
                        end
                        % Save data about the cell:
                        % Frame number
                        cellInfo{lane}(laneIndex(lane),1) = ii;
                        % Cell number
                        cellInfo{lane}(laneIndex(lane),2) = counter;
                        % Line intersection
                        cellInfo{lane}(laneIndex(lane),3) = line;
                        % Saves the area of the cell in pixels
                        cellInfo{lane}(laneIndex(lane),4) = cellCentroids(jj,1).Area(1);
                        % Store the length of the major axis
                        cellInfo{lane}(laneIndex(lane),5) = cellCentroids(jj,1).MajorAxisLength;
                        % Store the length of the minor axis
                        cellInfo{lane}(laneIndex(lane),6) = cellCentroids(jj,1).MinorAxisLength;
                        
                        % Updates the checking array and lane index
                        checkingArray(line,lane) = ii;
                        laneIndex(lane) = laneIndex(lane) + 1;
                        % Update workingFrame
                        workingFrame = workingFrame | currentRegion;
                    end
                end
            end
            % Sets line = 0 so if the cell is not on a line, it is not
            % counted next loop
            line = 0;
        end
    end
    
    %% Frame postprocessing
    % Save the labeled image
    % processedFrames(:,:,ii) = logical(workingFrame);
    
    if(WRITEVIDEO_FLAG)
        tempFrame = imoverlay(read(cellVideo,ii), bwperim(logical(workingFrame)), [1 1 0]);
        figure(2); imshow(tempFrame);
        writeVideo(outputVideo, tempFrame);
    end
    
    % Check to see if the arrays in 'cellInfo' are filling.  If there are
    % less than 10 more empty rows in any given array, estimate the number
    % of additional rows needed, based on the current filling and the
    % number of frames remaining.
    for jj = 1:16
        if((size(cellInfo{jj},2) - laneIndex(jj)) <= 10)
            vertcat(cellInfo{jj}, zeros(floor(((numFrames/ii-1)*size(cellInfo{jj},2))*1.1), 6));
        end
    end
    
%     % Progress bar update
%     if mod(ii, floor((numFrames)/100)) == 0
%         progressbar([],[], (ii/(numFrames)))
%     end
end

% Closes the video if it is open
if(WRITEVIDEO_FLAG)
    close(outputVideo);
end

%% Calls ProcessTrackingData to process the raw data and return
% transitData, an nx8 array where n is the number of cells that
% transited completely through the device.  The first column is the total
% transit time, the second gives the areas at each constriction, and columns
% 3-8 give the time taken to transit from constriction 1-2, 2-3, etc.

occl = zeros(size(processedFrames,3),17);
[occl, comp, debuglaneData, paired] = ProcessTrackingData(framerate, cellInfo, occl, currVideoName);
