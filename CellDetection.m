%% CellDetection.m
% function processed = CellDetection(currVideo, startFrame, endFrame, currPathName, currVideoName, mask, filters, lowContrast)
% CellDetection loads the videos selected earlier and processes each frame
% to isolate the cells.  When debugging, the processed frames can be
% written to a video, or the outlines of the detected cells can be overlaid
% on the video.

% Code from Dr. Amy Rowat's Lab, UCLA Department of Integrative Biology and
% Physiology
% Code originally by Bino Varghese (October 2011)
% Updated by David Hoelzle (January 2013)
% Updated by Mike Scott (July 2013)
% Rewritten by Ajay Gopinath (July 2013)
% Updated by Kendra Nyberg (May 2014)

% Inputs
%   - cellVideo: a videoReader object specifying a video to load
%   - startFrame: an integer specifying the frame to start analysis at
%   - endFrame: an integer specifying the frame to end analysis at
%   - folderName: a string specifying the filepath
%   - currVideoName: a string specifying the video's name
%   - mask: a logical array that was loaded in makeWaypoints and is used to
%       erase objects found outside of the lanes of the cell deformer.

% Outputs
%   - processed: An array of dimensions (height x width x frames) that
%       stores the processed frames.  Is of binary type.

% Changes
% Automation and efficiency changes made 03/11/2013 by Dave Hoelzle
% Commenting and minor edits on 6/25/13 by Mike Scott
% Increase in speed (~3 - 4x faster) + removed disk output unless debugging made on 7/5/13 by Ajay G.

function processed = CellDetection(currVideo, startFrame, endFrame, currVideoName, mask)

%%% This code analyzes a video of cells passing through constrictions
%%% to produce and return a binary array of the video's frames which
%%% have been processed to yield only the cells.

progressbar([],0,[])

DEBUG_FLAG = false; % flag for whether to show debug info
USEMASK_FLAG = true; % flag whether to binary AND the processed frames with the supplied mask
OVERLAYOUTLINE_FLAG = false; % flag whether to overlay detected outlines of cells on original frames
WRITEMOVIE_FLAG = false; % flag whether to write processed frames to a movie file

if(OVERLAYOUTLINE_FLAG)
    disp('!!Warning: OVERLAYOUTLINE_FLAG is set, frames cannnot be processed!!');
end


startTime1 = tic;

%% Initialization for debugging
% folderName = 'Y:\Kendra\HL60 Microfluidics\140803 - WT HL0 Drug Treatments I\WT\';
% currVideoName = 'dev5x5_200fps_2dppt_1o10PDMS_0,1F127_20x_0.4ms__35umfilter_001.avi';
% 
% currVideo = VideoReader([folderName, currVideoName]);
% startFrame = 1;
% endFrame = currVideo.NumberOfFrames;
% 
% % to determine mask variable
% [j,k] = regexp(currVideoName, 'dev\d*x'); % store start/end indices of template size
% templateSize = currVideoName((j+3):(k-1)); % removes 'dev' at the start, and 'x' at the end
% userPref = {'N','N','N', '1', 'N', '1'};
% % [mask, lineTemplate, xOffset, maskCheck] = MakeWaypoints(currVideo, templateSize, userPref(5), 'Y');
% [mask, lineTemplate, xOffset] = MakeWaypoints(currVideo, templateSize, userPref(5));

% if ismac
%     folderName = '/Volumes/Rowat Lab Data 2/Kendra/Lowry Lab Collaboration/150515 - Stem Cell Panel I/H1/';
%     currVideoName = 'dev5x5_200fps_24hppt_20x_0.4ms_030.avi';
%     
% else
%     folderName = 'Y:\Kendra\Lowry Lab Collaboration\150515 - Stem Cell Panel I\H1\';
%     currVideoName = 'dev5x5_200fps_24hppt_20x_0.4ms_030.avi';
% end
% 
% currVideo = VideoReader([folderName, currVideoName]);
% startFrame = 300;
% endFrame = currVideo.NumberOfFrames;
% 
% % Determines mask variable
% [j,k] = regexp(currVideoName, 'ev\d*x'); % store start/end indices of template size
% templateSize = currVideoName((j+2):(k-1)); % removes 'dev' at the start, and 'x' at the end
% userPref = {'N','N','N', '1', 'N', '1'};
% [mask, lineTemplate, xOffset, yOffset] = MakeWaypoints(currVideo, templateSize,  'N');
% 
% % Determines framerate
% [m, n] = regexp(currVideoName, '\d*fps'); % store start/end indices of frame rate
% frameRate = str2double(currVideoName(m:(n-3))); % removes 'fps'  at the end


%%
% Checks video format to ensure compatability
isVideoGrayscale = (strcmp(currVideo.VideoFormat, 'Grayscale') == 1);

disp(sprintf(['\nStarting cell detection for ', currVideoName, '...']));

% stores the number of frames that will be processed
effectiveFrameCount = (endFrame-startFrame+1) ;

% store the height/width of the cell video for clarity
height = currVideo.Height;
width = currVideo.Width;

% creates filters for cell detection
filters = [500; 3; 300; 25; 5];


%% Calculate initial background image
 sampleWindow = filters(1);

% if the sampling window is larger than the number of frames present,
% the number is set to all the frames present instead
if((sampleWindow+startFrame) > endFrame)
    sampleWindow = effectiveFrameCount-1;
end

% Store the first sampleWindow frames into bgFrames
% Store the first sampleWindow frames into bgFrames
bgFrames = zeros(height, width, sampleWindow, 'uint8');

if(isVideoGrayscale)
    bgFrames(:,:,startFrame:(startFrame+sampleWindow-1)) = read(currVideo, [startFrame (startFrame+sampleWindow-1)]); % store the frame that was read in bgFrames
else
    for j = startFrame:(startFrame+sampleWindow-1)
        temp = read(cellVideo, j);
        bgFrames(:,:,j) = temp(:,:,1);
    end
end


% calculate the initial 'background' frame for the first sampleWindow
% frames by storing the corresponding pixel value as the mean value of each
% corresponding pixel of the background frames in bgFrames

backgroundImg = uint8(mean(bgFrames, 3));
lastImg = read(currVideo, currVideo.NumberofFrames);
backgroundImg = imsubtract(backgroundImg, imsubtract(lastImg, backgroundImg));
lastBackgroundImg = double(backgroundImg);


% clear variables for better memory management
clear frameIdxs;



%% Prepare for Cell Detection
% create structuring elements used in cleanup of grayscale image
forClose = strel('disk', 5);

% automatic calculation of threshold value for conversion from grayscale to binary image
threshold = graythresh(backgroundImg) / 30;

% preallocate memory for marix for speed
if(OVERLAYOUTLINE_FLAG)
    processed = zeros(height, width, effectiveFrameCount, 'uint8');
else
    processed = false(height, width, effectiveFrameCount);
end

bgProcessTime = toc(startTime1);
startTime2 = tic;
%% Step through video
% iterates through each video frame in the range [startFrame, endFrame]
for frameIdx = startFrame:endFrame
    % reads in the movie file frame at frameIdx
    if(isVideoGrayscale)
        currFrame = read(currVideo, frameIdx);
    else
        temp = read(currVideo, frameIdx);
        currFrame = temp(:,:,1);
    end
    
    % if the current frame is after the first sampleWindow frames,
    % start adjusting the backgroundImage so that it represents a 'moving'
    % average of the pixel values of the frames in the interval
    % [frameIdx-sampleWindow, frameIdx]. This better localizes the background 
    % imageso it 'adapts' to the local frames and appears to better segment 
    % the cells. bgFrames is used to store the previous sampleWindow frames
    % so that memory is recycled.
    if(frameIdx >= sampleWindow+startFrame)
        bgImgDbl = lastBackgroundImg - double(bgFrames(:,:,(mod(frameIdx-1,sampleWindow)+1)))/sampleWindow + double(currFrame)/sampleWindow;
        backgroundImg = uint8(bgImgDbl);
        backgroundImg = imsubtract(backgroundImg, imsubtract(lastImg, backgroundImg));
        lastBackgroundImg = bgImgDbl;
        bgFrames(:,:,mod(frameIdx-1,sampleWindow)+1) = currFrame;
    end
    
    
    if (frameIdx == 1)
        % Do cell detection for first frame to determine number of
        % occlusions
        % converts first frame to bw with threshold + without background subtraction
        cleanImg = im2bw(currFrame(1:60,:), graythresh(currFrame)*1.1);
        cleanImg = imbothat(cleanImg, forClose);
        cleanImg = bwareaopen(cleanImg, 35);
        cleanImg = imclose(cleanImg, forClose);
        cleanImg = imfill(cleanImg, 'holes')-cleanImg;
        cleanImg = vertcat(cleanImg, zeros(size(currFrame,1)-60, size(currFrame,2)));

    else
        % Do cell detection
        % subtracts the background in bgImgs from each frame, hopefully leaving
        % just the cells
        cleanImg = imsubtract(backgroundImg, currFrame);
        cleanImg = im2bw(cleanImg, threshold); % thresholds
        
        % Cleanup 
        % clean the grayscale image of the cells to improve detection
        cleanImg = bwareaopen(cleanImg, 20); %removes small particles
        cleanImg = medfilt2(cleanImg, [2, 2]); %median smoothing
        cleanImg = imclose(cleanImg, forClose); %connects gaps
        cleanImg = imfill(cleanImg, 'holes'); %fills in holes
        cleanImg = bwareaopen(cleanImg, 130); %removes small particles
        
%         figure(2); imshow(horzcat(currFrame, cleanImg*255));
    end 
   
    
    if(USEMASK_FLAG)
        % binary 'AND' to find the intersection of the cleaned image and the mask
        % to prevent the detected cell boundaries from being outside the microfluidic device
        cleanImg = cleanImg & mask; 
    end
    
    if OVERLAYOUTLINE_FLAG == 1
        cleanImg = imadd(currFrame, uint8(bwperim(cleanImg)*255));
%         cleanImg = imoverlay(cleanImg, logical(lineTemplate),[1 1 0]);
    end
    
    %% Store cleaned image of segmented cells in processed
    processed(:,:,frameIdx-startFrame+1) = cleanImg;
 
    
    % Increments the progress bar, each time 1% of the frames are finished
    if mod(frameIdx, floor(effectiveFrameCount/100)) == 0
        progressbar([], frameIdx/effectiveFrameCount, [])
    end
end

% stop recording the time and output debugging information
framesTime = toc(startTime2);
totalTime = bgProcessTime + framesTime;
disp(['Time taken for cell detection: ', num2str(totalTime), ' secs']);
disp(['Average time to detect cells per frame: ', num2str(framesTime/effectiveFrameCount), ' secs']);

%% Set up frame viewer and write to file if debugging is on
if(DEBUG_FLAG)
    implay(processed);
    
    % if video file is set
    if(WRITEMOVIE_FLAG)
        writer = VideoWriter([folderName, 'cellsdetected_', videoName]);
        open(writer);
        
        if(islogical(processed))
            processed = uint8(processed); % convert to uint8 for use with writeVideo

            % make binary '1's into '255's so all resulting pixels will be
            % either black or white
            for idx = 1:effectiveFrameCount
                processed(:,:,idx) = processed(:,:,idx)*255;
            end
        end
        
        % write processed frames to disk
        for currFrame = startFrame:endFrame
            writeVideo(writer, processed(:,:,currFrame-startFrame+1));
        end
        
        close(writer);
    end
end
