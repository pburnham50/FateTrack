classdef FTselector < handle

    properties (Access = public)

        pointTableHandle % points in Table
        imageToFrameTable
        fileTable % all the files and their frameNumbers

        currPoints % points in image

        %imageToFrameTable % This should link the image filenames to frames
        imageFiles % Can probably fold this into the above.



        figHandle

        pointController

        lassoButtonHandle
        deselectAllButtonHandle
        connectParentButtonHandle
        saveButtonHandle
        showPointIDsHandle
        showMasksHandle
        showCytoMasksHandle
        currentFramePopupHandle

        noneButtonHandle
        cyfivetwoButtonHandle
        cyfiveoneButtonHandle
        cythreetwoButtonHandle
        yfptwoButtonHandle
        yfponeButtonHandle
        yfponecyfivetwoButtonHandle
        latentButtonHandle

        axesHandle
        imageHandle
        linesHandle

    end

    methods
        function p = FTselector(inFilename)

            p.figHandle = figure('Visible','off','Position',[360,500,450,285]);
            p.figHandle.KeyPressFcn = @p.GUIWindowKeyPressFcn;

            p.lassoButtonHandle = uicontrol('Style','pushbutton','String','LassoSelect','Position',[5,250,80,25]);
            p.deselectAllButtonHandle = uicontrol('Style','pushbutton','String','deselect','Position',[5,280,80,25]);
            p.connectParentButtonHandle = uicontrol('Style','pushbutton','String','region','Position',[5,310,80,25]);
            p.saveButtonHandle = uicontrol('Style','pushbutton','String','save','Position',[5,340,80,25]);
            p.showPointIDsHandle = uicontrol('Style','checkbox','String','showID','Position',[5,370,80,25]);
            p.showMasksHandle = uicontrol('Style','checkbox','String','showMask','Position',[5,400,80,25]);
            p.showCytoMasksHandle = uicontrol('Style','checkbox','String','showCytoMask','Position',[5,430,80,25]);

            p.noneButtonHandle =        uicontrol('Style','pushbutton','String','none','Position',[5,460,80,25]);
            p.yfptwoButtonHandle =         uicontrol('Style','pushbutton','String','2-yfp','Position',[5,490,80,25]);
            p.cyfivetwoButtonHandle =        uicontrol('Style','pushbutton','String','2-cy5','Position',[5,520,80,25]);
            p.cythreetwoButtonHandle =   uicontrol('Style','pushbutton','String','2-cy3','Position',[5,550,80,25]);
            p.yfponeButtonHandle =      uicontrol('Style','pushbutton','String','1-yfp','Position',[5,580,80,25]);
            p.cyfiveoneButtonHandle =     uicontrol('Style','pushbutton','String','1-cy5','Position',[5,610,80,25]);
            p.yfponecyfivetwoButtonHandle = uicontrol('Style','pushbutton','String','1-yfp_2-cy5','Position',[5,640,80,25]);
            p.latentButtonHandle =   uicontrol('Style','pushbutton','String','latent','Position',[5,670,80,25]);

            p.currentFramePopupHandle = uicontrol(p.figHandle,'Style','popupmenu');
            p.currentFramePopupHandle.Position = [5 5 200 25];




            % Use following
            %currFrame = str2num(popupHandle.String{popupHandle.Value});

            p.axesHandle = axes('Units','normalized','Position',[0.05 0.02 .98 .98]);

            imshow(rand(100),'Parent',p.axesHandle);

            % Maybe keep invisible until we load data?
            p.figHandle.Visible = 'on';

            % Should fix this with a proper display method

            p.pointController = pointController(p.pointTableHandle, p.imageToFrameTable, p.deselectAllButtonHandle, p.lassoButtonHandle, p.connectParentButtonHandle, p.currentFramePopupHandle, p.axesHandle, p.linesHandle);
            pc = p.pointController;

            p.lassoButtonHandle.Callback ={@pc.lassoButtonPushed}; % Needs fixing
            p.deselectAllButtonHandle.Callback ={@pc.deselectAllButtonPushed}; % Needs fixing
            p.connectParentButtonHandle.Callback ={@pc.connectParentButtonPushed}; % Needs fixing
            p.currentFramePopupHandle.Callback ={@pc.updateFrame}; % Needs fixing
            p.saveButtonHandle.Callback ={@pc.saveButtonPushed}; % Needs fixing
            p.showPointIDsHandle.Callback ={@pc.showPointIDsPushed}; % Needs fixing
            p.showMasksHandle.Callback ={@pc.showMasksPushed}; % Needs fixing
            p.showCytoMasksHandle.Callback ={@pc.showCytoMasksPushed}; % Needs fixing

            p.noneButtonHandle.Callback ={@pc.setAnnotation}; % Needs fixing
            p.cyfivetwoButtonHandle.Callback ={@pc.setAnnotation}; % Needs fixing
            p.cyfiveoneButtonHandle.Callback ={@pc.setAnnotation}; % Needs fixing
            p.cythreetwoButtonHandle.Callback ={@pc.setAnnotation}; % Needs fixing
            p.yfptwoButtonHandle.Callback ={@pc.setAnnotation}; % Needs fixing
            p.yfponeButtonHandle.Callback ={@pc.setAnnotation}; % Needs fixing
            p.latentButtonHandle.Callback ={@pc.setAnnotation}; % Needs fixing
            p.yfponecyfivetwoButtonHandle.Callback ={@pc.setAnnotation};

            pc.showPointIDsHandle = p.showPointIDsHandle;
            pc.showMasksHandle = p.showMasksHandle;
            pc.showCytoMasksHandle = p.showCytoMasksHandle;
            pc.saveFilename = inFilename;


            p.loadAllData();
            pc.pointTableHandle = p.pointTableHandle;
            pc.imageToFrameTable = p.imageToFrameTable;
            %p.pointController.saveButtonPushed([],[]);

            pc.showImages();
            pc.makePoints();

        end


        function p = loadAllData(p) % Lots of this should probably be in a fileController

            %%% SHOULD USE THIS
            p.fileTable = parseFiles();
            %Probably should call writetable within parseFiles
            writetable(p.fileTable,'fileTable.csv');

            % Replace with some sort of file loader
            p.imageFiles = p.fileTable.fileName;

            p.currentFramePopupHandle.Value = 1;
            p.currentFramePopupHandle.String = p.imageFiles(1:end-1);
            p.currentFramePopupHandle.UserData = p.imageFiles;

            p.pointTableHandle = pointTable(p.pointController.saveFilename);

        end


        function p = GUIWindowKeyPressFcn(p, src, eventdata)
            % determine the key that was pressed
            keyPressed = eventdata.Key;
            switch(keyPressed)
                case 'd'
                    p.pointController.deselectAllButtonPushed(src, eventdata);
                case 'f'
                    p.pointController.deleteButtonPushed(src, eventdata);
                case 'g'
                    p.pointController.connectParentButtonPushed(src, eventdata);
                case 'c'
                    p.currentFramePopupHandle.Value = max([1 p.currentFramePopupHandle.Value-1]);
                    p.pointController.updateFrame(src, eventdata);
                case 'v'
                    p.currentFramePopupHandle.Value = min([length(p.currentFramePopupHandle.String) p.currentFramePopupHandle.Value+1]);
                    p.pointController.updateFrame(src, eventdata);
                case 'z'
                    p.pointController.zoomMode();
                case 'q'
                    p.pointController.zoomReset();
                case 'x'
                    p.pointController.unZoom();
                case 't'
                    p.pointController.showPointIDsHandle.Value = ~p.pointController.showPointIDsHandle.Value;
                    p.pointController.showPointIDsPushed(p.pointController.showPointIDsHandle,eventdata);
            end

        end

    end
end
