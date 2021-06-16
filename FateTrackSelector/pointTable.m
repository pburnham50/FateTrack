classdef pointTable < handle

    properties (Access = public)

        allPoints

    end

    methods

        function p = pointTable(varargin)
            if nargin == 0
                fprintf('Please Add Table\n');
            else
                fprintf('Loading Table\n');
                p.allPoints = readtable(varargin{1},'TextType','string');
                if ~ismember("annotation",p.allPoints.Properties.VariableNames)
                   col = repmat("none",height(p.allPoints),1);
                   p.allPoints.annotation = col;
                end
            end
        end

        function p = changeAnnotation(p, pointIDs, newAnnotation)
            idx = ismember(p.allPoints.pointID,pointIDs);
            p.allPoints.annotation(idx) = newAnnotation;
        end

        % remove add raw points

        %remove removePoints

        %remove guessParent

        %remove guessParents

        %remove assignParents


        function p = setPointCoordinates(p, pointID, xCoord, yCoord)
            p.allPoints(p.allPoints.pointID == pointID,:).xCoord = xCoord;
            p.allPoints(p.allPoints.pointID == pointID,:).yCoord = yCoord;
        end

        %remove assignParents

        function outTab = getAllPointsInFrame(p)
            % outTab = p.allPoints(p.allPoints.frameNumber == frameNumber,:);
            outTab = p.allPoints;
        end

        % Probably need something to set an entire frame's worth of points
        % UNTESTED!! This also might be done better in a more abstract way
        % Also, not clear this is actually needed.'
        function p = updateEntireFrame(p,frameNumber,newTable)
            idx = p.allPoints.frameNumber == frameNumber;
            p.allPoints(idx,:) = [];
            p.allPoints = [p.allPoints newTable];
        end


    end

end
