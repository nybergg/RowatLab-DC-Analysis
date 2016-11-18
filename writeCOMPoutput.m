%% writeSPSSoutput.m
%This function formats data in a way that is easily transferable to the
%statistics program SPSS.

% Code from Dr. Amy Rowat's Lab, UCLA Department of Integrative Biology and
% Physiology
% Code by Sam Bruce (July 2013)

function [writeData, folderNames] = writeCOMPoutput(pathNames, writeData, i,folderNames, compData)

%Analazes every selected pathname for their unique features. Deletes all
%common foldernames except for one.
%Initializes writeToSPSS with its largest horizontal size for selected
%videos.

if (isempty(folderNames) && isempty(writeData))
    for ii = 1:size(pathNames,2)
        pathlength(ii)= size(strfind(pathNames{ii},'\'),2);
    end
    
    folderNames = cell(0,max(pathlength));
    for ii = 1:size(pathNames,2)
        pathName = pathNames{ii};
        folderName = strsplit(pathName(1:size(pathName,2)-1),'\\');
        folderNames(ii,1:size(folderName,2)) = folderName;
    end
    
    
    while (size(folderNames,2)>1 & strcmp(folderNames{1,2},folderNames(:,2)))
        folderNames(:,1) = [];
    end
    
    writeData = cell(0,9+size(folderNames,2));
    
end

%Adds each folder name in the file's path to the cell array for future
%sorting in SPSS.
folderName = folderNames(i,:);

for jj = 1:size(compData,1)
    for k=1:size(folderName,2)
         compData{jj,10+k} = folderName{k};
    end
end

%Vertically concatinates new data into writeToSPSS
if(~isempty(writeData))
    writeData = vertcat(writeData, compData);
else
    writeData = compData;
end

