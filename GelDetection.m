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

function processed = GelDetection(currVideo, startFrame, endFrame, currVideoName, mask, filters)

%%% This code analyzes a video of cells passing through constrictions
%%% to produce and return a binary array of the video's frames which
%%% have been processed to yield only the cells.

progressbar([],0,[])


DEBUG_FLAG = false; % flag for whether to show debug info
WRITEMOVIE_FLAG = false; % flag for whether to write processed frames to movie on disk
USEMASK_FLAG = true; % flag whether to binary AND the processed frames with the supplied mask
OVERLAYOUTLINE_FLAG = false; % flag whether to overlay detected outlines of cells on original frames

if(OVERLAYOUTLINE_FLAG)
    disp('!!Warning: OVERLAYOUTLINE_FLAG is set, frames cannnot be processed!!');
end


startTime1 = tic;

%% Initialization for debugging
folderName = 'Y:\Kendra\Agarose Microgels\160504 - Agarose Microgels with Varying Concentrations Batch II - I\1.0\';
currVideoName = 'dev5x5_200fps_4psi_0,1percF127_20x_000.avi';

currVideo = VideoReader([folderName, currVideoName]);
startFrame = 1;
endFrame = currVideo.NumberOfFrames;

% to determine mask variable
[j,k] = regexp(currVideoName, 'dev\d*x'); % store start/end indices of template size
templateSize = currVideoName((j+3):(k-1)); % removes 'dev' at the start, and 'x' at the end
userPref = {'N','N','N', '1', 'N'};
% [mask, lineTemplate, xOffset, maskCheck] = MakeWaypoints(currVideo, templateSize, userPref(5), 'Y');
[mask, lineTemplate, xOffset] = MakeWaypoints(currVideo, templateSize, userPref(5));

filters = [500; 3; 140; 35; 5];


%%


isVideoGrayscale = (strcmp(currVideo.VideoFormat, 'Grayscale') == 1);

disp(sprintf(['\nStarting cell detection for ', currVideoName, '...']));

% stores the number of frames that will be processed
effectiveFrameCount = (endFrame-startFrame+1) ;

% store the height/width of the cell video for clarity
height = currVideo.Height;
width = currVideo.Width;



%% Calculate initial background image
 sampleWindow = filters(1);


% if the sampling window is larger than the number of frames present,
% the number is set to all the frames present instead
if((sampleWindow+startFrame) > endFrame)
    sampleWindow = effectiveFrameCount-1;
end

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


for ii = 1:size(bgFrames, 3)
    bgFrames(:,:,ii) =  imadjust(bgFrames(:,:,ii),[0.32; 1.00],[0.00; 1.00], 0.30);
end


% calculate the initial 'background' frame for the first sampleWindow
% frames by storing the corresponding pixel value as the mean value of each
% corresponding pixel of the background frames in bgFrames
backgroundImg = uint8(mean(bgFrames, filters(2)));

% clear variables for better memory management
clear frameIdxs;

%% Prepare for Cell Detection
% create structuring elements used in cleanup of grayscale image
forClose = strel('disk', 10);

% automatic calculation of threshold value for conversion from grayscale to binary image
threshold = graythresh(backgroundImg) / filters(3);


% preallocate memory for marix for speed
if(OVERLAYOUTLINE_FLAG)
    processed = zeros(height, width, effectiveFrameCount, 'uint8');
else
    processed = false(height, width, effectiveFrameCount);
end


bgProcessTime = toc(startTime1);
startTime2 = tic;
lastBackgroundImg = double(backgroundImg);


%% Gel detection
% iterates through each video frame in the range [startFrame, endFrame]
for frameIdx = startFrame:2:endFrame
    % reads in the movie file frame at frameIdx
    if(isVideoGrayscale)
        currFrame = imadjust(read(currVideo, frameIdx),[0.32;
            1.00],[0.00; 1.00], 0.30); %Agarose Presettings
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
        lastBackgroundImg = bgImgDbl;
        bgFrames(:,:,mod(frameIdx-1,sampleWindow)+1) = currFrame;
    end
    
    % Gel detection
%     cleanImg = im2bw(medfilt2(imsubtract(backgroundImg, currFrame), [3, 3]), threshold);
%     figure(2); imshow(cleanImg*255)
    cleanImg = im2bw(imsubtract(backgroundImg, currFrame), threshold);


    % Cleanup
    % clean the grayscale image of the gels to improve detection
    cleanImg = medfilt2(cleanImg, [filters(5), filters(5)]);
    cleanImg = bwareaopen(cleanImg, filters(4));
    cleanImg = imfill(cleanImg, 'holes');
    cleanImg = imtophat(cleanImg, forClose);
    cleanImg = imclose(cleanImg, forClose);
    

    
    if(USEMASK_FLAG)
        % binary 'AND' to find the intersection of the cleaned image and the mask
        % to prevent the detected cell boundaries from being outside the microfluidic device
        cleanImg = cleanImg & mask;
    end
    
    
    
    if OVERLAYOUTLINE_FLAG == 1
        cleanImg = imadd(currFrame, uint8(bwperim(cleanImg)*255));
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
        writer = VideoWriter([currPathName, 'cellsdetected_', currVideoName]);
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
