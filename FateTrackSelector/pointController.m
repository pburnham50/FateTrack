classdef pointController < handle

    properties (Access = public)

        pointTableHandle % links points in Table
        currPoints % links points in image to pointIDs

        zoomObject % for switching zoom on and off.

        imageToFrameTable % This should link the image filenames to frames

        saveFilename

        GFPStatus = false;
        TransStatus = false;

        lassoButtonHandle
        deselectAllButtonHandle
        connectParentButtonHandle
        currentFramePopupHandle
        showPointIDsHandle
        showMasksHandle
        showPointIDsWithInfoHandle
        axesHandle
        linesHandle
        imageHandle

    end

    methods

        % Constructor will "register" components
        % This could be more elegant, because needs to be updated every
        % time. Probably could do something with events and listeners,
        % whatev.
        function p = pointController(pointTableHandle, imageToFrameTable, deselectAllButtonHandle, lassoButtonHandle, connectParentButtonHandle, currentFramePopupHandle, axesHandle, linesHandle)
            p.pointTableHandle = pointTableHandle;
            p.imageToFrameTable = imageToFrameTable;

            p.lassoButtonHandle = lassoButtonHandle;
            p.deselectAllButtonHandle = deselectAllButtonHandle;
            p.connectParentButtonHandle = connectParentButtonHandle;
            p.currentFramePopupHandle = currentFramePopupHandle;
            p.axesHandle = axesHandle;
            p.linesHandle = linesHandle;
        end

        function p = zoomMode(p,src,eventdata)
            %p.zoomObject = zoom;
            %p.zoomObject.Enable = 'on';
            zoom(1.5);


            % The following enables you to add key functions, but thats
            % not what we want.
            %hManager = uigetmodemanager(gcf);
            %[hManager.WindowListenerHandles.Enabled] = deal(false);

            %figHandle = gcf;
            %figHandle.KeyPressFcn = @tempWindowKeyPressFcn;
            %p.zoomObject.ActionPostCallback = @p.undoZoom;
        end

        function p = unZoom(p,src,eventdata)
            zoom(1/1.5);
        end

        function p = zoomReset(p,src,eventdata)
            zoom out;
        end

        % unused
        function p = undoZoom(p,src,eventdata)
            hManager = uigetmodemanager(gcf);
            [hManager.WindowListenerHandles.Enabled] = deal(false);
            p.zoomObject.Enable = 'off';
        end

        % This makes points
        function p = makePoints(p) % Probably should be in the view but whatever.
          %  popupHandle = p.currentFramePopupHandle;
          %  currFrame = popupHandle.Value; %str2num(popupHandle.String{popupHandle.Value});
            Tcurr = p.pointTableHandle.allPoints ;

            if ~isempty(p.currPoints)
                if isvalid(p.currPoints)
                    delete(p.currPoints);
                end
            end


            p.currPoints = images.roi.Point;
            for i = 1:height(Tcurr)
                % make point a different color if it is already annotated
                % check for annotation
                idx = p.pointTableHandle.allPoints.pointID == Tcurr.pointID(i);
                annotation = p.pointTableHandle.allPoints.annotation(idx);
                if annotation ~= "none" % different color for annotated point
                    p.currPoints(i) = drawpoint(p.axesHandle,'Position',[Tcurr.xCoord(i) Tcurr.yCoord(i)],...
                        'Color','magenta','SelectedColor','c');
                else
                    p.currPoints(i) = drawpoint(p.axesHandle,'Position',[Tcurr.xCoord(i) Tcurr.yCoord(i)],...
                        'Color',[.4 .4 1],'SelectedColor','c');
                end
                p.currPoints(i).UserData = Tcurr.pointID(i);
                p.currPoints(i).Label = getPtLabel(p,Tcurr.pointID(i));
            end

        end

        % removed drawLines

        % removed addCurrPointButtonPushed

        % removed addNextPointButtonPushed

        % removed addNextConnectPointButtonPushed

        function p = saveButtonPushed(p,src,eventdata)
            fprintf('SAVE\n');
            writetable(p.pointTableHandle.allPoints,p.saveFilename);
            % p.pointTableHandle.allPoints.writetable(p.saveFilename);
        end

        function p = setAnnotation(p,src,eventdata)
            currPtSelected = find([p.currPoints(:).Selected]);

            pointIDs = [p.currPoints(currPtSelected).UserData];

            annotation = src.String;

            p.pointTableHandle.changeAnnotation(pointIDs,annotation);
            p.showPointIDsPushed(src,eventdata);
            p.deselectAllButtonPushed(src,eventdata);


        end

        % Probably could make this abstract or whatever kind of method
        function label = getPtLabel(p,theID)
            %%% NEEDS UPDATING TO PULL THE LABEL FROM THE POINTTABLE
            if p.showPointIDsHandle.Value
                idx = p.pointTableHandle.allPoints.pointID == theID;
                annotation = p.pointTableHandle.allPoints.annotation(idx);
                if annotation ~= "none"
                    label = num2str(theID)+" "+annotation;
                else
                    label = num2str(theID);
                end
            elseif p.showPointIDsWithInfoHandle.Value
                idx = p.pointTableHandle.allPoints.pointID == theID;
                annotation = p.pointTableHandle.allPoints.annotation(idx);
                channel_file = char(p.currentFramePopupHandle.UserData{p.currentFramePopupHandle.Value});
                channel_string = channel_file(1:end-4);
                lp_column = "lp_Rd" + channel_string(1) + "_" + channel_string(2:end) + "_mean_intensity";
                nuc_column = "nuc_Rd" + channel_string(1) + "_" + channel_string(2:end) + "_mean_intensity";
                hcr_val_str = "LPe5=" + p.pointTableHandle.allPoints.(lp_column)(idx) + ", nuc=" + p.pointTableHandle.allPoints.(nuc_column)(idx);
                if annotation ~= "none"
                    label = num2str(theID)+" "+annotation + " " + hcr_val_str;
                else
                    label = num2str(theID) + " " + hcr_val_str;
                end
            else
                label = '';
            end
        end

        function p = getMask(p)
            %%% NEEDS UPDATING TO PULL THE LABEL FROM THE POINTTABLE
            if p.showMasksHandle.Value
                p.showImagesNoPointUpdateWithMask();
            else
                p.showImagesNoPointUpdate();
            end
        end
        
%         function p = getCytoMask(p)
%             %%% NEEDS UPDATING TO PULL THE LABEL FROM THE POINTTABLE
%             if p.showCytoMasksHandle.Value
%                 p.showImagesNoPointUpdateWithCytoMask();
%             else
%                 p.showImagesNoPointUpdate();
%             end
%         end

        function p = showPointIDsPushed(p,src,eventdata)
            for i = 1:length(p.currPoints)
                if ~isempty(p.currPoints(i).UserData)
                    p.currPoints(i).Label = getPtLabel(p,p.currPoints(i).UserData);
                end
            end
        end
        
        function p = showPointIDsWithInfoPushed(p,src,eventdata)
            for i = 1:length(p.currPoints)
                if ~isempty(p.currPoints(i).UserData)
                    p.currPoints(i).Label = getPtLabel(p,p.currPoints(i).UserData);
                end
            end
        end

        function p = showMasksPushed(p,src,eventdata)
            p.getMask()
        end
        
%         function p = showCytoMasksPushed(p,src,eventdata)
%             p.getCytoMask()
%         end

        function p = deselectAllButtonPushed(p,src,eventdata)
            for i = 1:length(p.currPoints)
                p.currPoints(i).Selected = 0;
            end
        end

        % removed connectParentButtonPushed

        function p = lassoButtonPushed(p,src,eventdata)
            % Basically, send messages to deletePoint for each selected. Then, clean
            % up the list.
            % Then, drawLines

            for i = 1:length(p.currPoints)
                pt = p.currPoints(i);
                if pt.Selected
                    p.deletePoint(pt,[]);
                end
            end

            p.currPoints = p.currPoints(isvalid(p.currPoints));
            p.showPointIDsPushed(src,eventdata);
            p.drawLines();

        end

        %removed deletePoint(p,src,eventdata)

        %removed pointMoved

        function p = updateFrame(p,src,eventData)
            p.showImages();
            fprintf('SAVE\n');
            writetable(p.pointTableHandle.allPoints,p.saveFilename);
        end

        %removed toggleGFP

        %removed toggleTrans

        function p = showImagesNoPointUpdate(p)
            currFrame = p.currentFramePopupHandle.Value;
            p.imageToFrameTable = parseFiles()
            image1 = imread(p.currentFramePopupHandle.UserData{currFrame});
            image1 = imadjust(image1,[],[],0.1);
            toggleImage = zeros([size(image1),3]);

            outRGB = makeColoredImage(scale(im2double(image1)),[0.9648 0.5781 0.1172]) + toggleImage;

            if ~isempty(p.imageHandle)
                p.imageHandle.CData = outRGB;
            else
                p.imageHandle = imshow(outRGB,'Parent',p.axesHandle);
            end
        end

        function p = showImagesNoPointUpdateWithMask(p)
            currFrame = p.currentFramePopupHandle.Value;
            p.imageToFrameTable = parseFiles()
            image1 = imread(p.currentFramePopupHandle.UserData{currFrame});
            image1 = imadjust(image1,[],[],0.1);
            image2 = imread(p.imageToFrameTable(p.imageToFrameTable.wavelength=="mask",:).fileName{end});
            toggleImage = zeros([size(image1),3]);

            outRGB = makeColoredImage(scale(im2double(image2)),[0 0.6797 0.9336]) + makeColoredImage(scale(im2double(image1)),[0.9648 0.5781 0.1172]) + toggleImage;

            if ~isempty(p.imageHandle)
                p.imageHandle.CData = outRGB;
            else
                p.imageHandle = imshow(outRGB,'Parent',p.axesHandle);
            end

        end
        
%         function p = showImagesNoPointUpdateWithCytoMask(p)
%             currFrame = p.currentFramePopupHandle.Value;
%             p.imageToFrameTable = parseFiles()
%             image1 = imread(p.currentFramePopupHandle.UserData{currFrame});
%             image1 = imadjust(image1,[],[],0.1);
%             image2 = imread(p.imageToFrameTable(p.imageToFrameTable.wavelength=="cyto_mask",:).fileName{end});
%             toggleImage = zeros([size(image1),3]);
% 
%             outRGB = makeColoredImage(scale(im2double(image2)),[0 0.6797 0.9336]) + makeColoredImage(scale(im2double(image1)),[0.9648 0.5781 0.1172]) + toggleImage;
% 
%             if ~isempty(p.imageHandle)
%                 p.imageHandle.CData = outRGB;
%             else
%                 p.imageHandle = imshow(outRGB,'Parent',p.axesHandle);
%             end
% 
%         end

        function p = showImages(p)

            % p.makePoints();
            p.showImagesNoPointUpdate();
            % bringToFront(p.currPoints)
        end

        function p = showPoints(p)

            p.makePoints();

        end


        % Make a button to deselect all.

        % make a method to update the connection between the lines as we
        % move the child/parent. Need to make this only affect the relevant
        % line object so we don't have to draw them all again.

        % This is how to use a method for a callback. obj is the instance.
        % uicontrol('Style','slider','Callback',@obj.sliderCallback)
        % function methodName(obj,src,eventData)
    end
end
