
function fileTable = parseFiles()

dd = dir('*.tif');

fileTable = cell2table(cell(0,3), 'VariableNames', {'fileName','round','wavelength'});

for i = 1:length(dd)
    currFile = dd(i).name;
    [filepath,name,ext] = fileparts(currFile);
    C = strsplit(name,'_');
    if C{end}~="mask"
      roundnum = C{end}(1);
      wl = C{end}(2:end);
    else
        if C{1} == "cyto"
            roundnum = "cyto_mask";
            wl = "cyto_mask";
        else
            roundnum = "mask";
            wl = "mask";
        end
    end

    fileTable = [fileTable;{currFile,roundnum,wl}];

end

%uniqueTimes = sort(unique(fileTable.time));

%frameTable = table(uniqueTimes,(1:length(uniqueTimes))','VariableNames',{'time','frameNumber'});


%fileTable = join(fileTable,frameTable);

%idx = fileTable.wavelength == "cy5";
