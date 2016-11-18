function [laneData] = CleanData(laneData)

% Eliminates ghost cells from background subtractions
    if (~isempty(laneData))
        % Finds rows in laneData where cells transit all the way to the 8th
        % line
        [compTransRow, ~] = find(laneData(:,7,1) ~= 0);
        % Finds rows in laneData where cells are only detected in the first
        % constriction
        [c1Row, ~] = find(laneData(:,4,1) == 0);

        
        % Excludes ghost images
        
        % Finds ghost images from background detection: have to meet
        % three conditions: 1) arrives at the first constriction due to
        % a real, fully transiting cell (i.e. the row after a
        % compTransRow row. 2) arrives before the real cell reaches the
        % eighth lane 3) never leaves the first constriction.
                  
        if (~isempty(c1Row) && ~isempty(compTransRow))
            if (~isequal(c1Row(1),1))
                [xStart, ~] = find(laneData(c1Row,2,1) < laneData(c1Row-1,7,1));
            else
                [xStart, ~] = find(laneData(c1Row(2:end),2,1) < laneData(c1Row(2:end)-1,7,1));
                
                if ~isempty(xStart)
                    xStart = xStart + 1;
                end
                
            end
            
            
            if (~isempty(xStart))
                
                xStart = c1Row(xStart);
                for b = 1:size(xStart, 1)
                    [xEnd, ~]  = find(compTransRow > xStart(b), 1, 'first');
                    if (~isempty(xEnd))
                        xEnd = compTransRow(xEnd);                        
                        
                        % Eliminate ghost cells
                        laneData = vertcat(laneData(1:xStart(b)-1, :, :), laneData(xEnd:end, :, :));
                        
                        % Eliminate indexing of ghost cells
                        % laneData row index of c1 cell (elimination and
                        % adjustment)
                        c1Row = vertcat(c1Row(c1Row < xStart(b)), c1Row(c1Row >= xEnd) - (xEnd - xStart(b)));
                        
                        % laneData row index of completely transiting cell
                        % (elimination and adjustment
                        compTransRow = vertcat(compTransRow(compTransRow < xEnd), compTransRow(compTransRow > xEnd) - (xEnd - xStart(b)));
                        
                        % Adjust consequetive xStart values
                        xStart =  xStart - (xEnd - xStart(b));
                        
                    else
                        % Eliminate the last laneData rows if they
                        % correspond to a ghost image.
                        laneData = laneData(1:xStart(b)-1, :, :);
                        
                        % Eliminate indexing of ghost cells
                        % laneData row index of c1 cell (elimination and
                        % adjustment)
                        c1Row = c1Row(c1Row < xStart(b));        
                        
                    end
                end
            end
        end
        
        

        
       % Reconnects lost cells when they dissappear in the detection when
       % they stay too long in the first constriction.

       if (~isempty(c1Row)) % checks for cells stuck in C1
            if (~isequal(c1Row(1),1)) 
                [xStart, ~] = find(laneData(c1Row-1,7,1) > 0); % finds cells that enter after fully transiting cells to determine the first detected frame 

            else
                % When the cell is the first cell to passage through
                [xStart, ~] = find(laneData(c1Row(2:end)-1,7,1) > 0);
                
                if ~isempty(xStart)
                    xStart = vertcat(1, xStart+1); %Adds the first cell if needed
                else
                    xStart = 1; % if the first cell is the only cell
                end

            end
            
            if (~isempty(xStart)) % if stuck cells have been detected
                
                xStart = c1Row(xStart); % identifies the beginning of cells stuck in C1
                
                if (~isempty(compTransRow)) 
                    
                    for b = 1:size(xStart, 1)
                        
                        [xEnd, ~]  = find(compTransRow > xStart(b), 1, 'first');
                        
                        if (~isempty(xEnd))
                            xEnd = compTransRow(xEnd);
                            
                            if ~isequal(xStart(b)+1,xEnd)
                                laneData(xEnd, 1:2, :) = laneData(xStart(b), 1:2, :); % setting the frame number at line 1 and 2 as correct start frames
                            end
                            
                            if xStart(b) ~= 1
                                % Eliminates rows in which the cell was
                                % stuck
                                laneData = vertcat(laneData(1:xStart(b)-1, :, :), laneData(xEnd:end, :, :));

                                % Eliminate indexing of ghost cells
                                % laneData row index of c1 cell (elimination and
                                % adjustment)
                                c1Row = vertcat(c1Row(c1Row < xStart(b)), c1Row(c1Row >= xEnd) - (xEnd - xStart(b)));

                                % laneData row index of completely transiting cell
                                % (elimination and adjustment
                                compTransRow = vertcat(compTransRow(compTransRow < xEnd), compTransRow(compTransRow >= xEnd) - (xEnd - xStart(b)));

                                % Adjust consequetive xStart values
                                xStart = xStart - (xEnd - xStart(b));

                                
                            else
                                laneData = laneData(xEnd:end, :, :);
                                
                                % Eliminate indexing of ghost cells
                                % laneData row index of c1 cell (elimination and
                                % adjustment)
                                c1Row = c1Row(c1Row >= xEnd) - (xEnd - xStart(b));

                                % laneData row index of completely transiting cell
                                % (elimination and adjustment
                                compTransRow = compTransRow(compTransRow >= xEnd) - (xEnd - xStart(b));

                                % Adjust consequetive xStart values
                                xStart = xStart - (xEnd - xStart(b));
                                
                            end 

                        else
                            laneData = laneData(1:xStart(b), :, :);
                            
                            % Eliminate indexing of ghost cells
                            % laneData row index of c1 cell (elimination and
                            % adjustment)
                            break
                        end
                    end
                else
                    laneData = laneData(xStart,:,:);
                end
            end            
       end
    end