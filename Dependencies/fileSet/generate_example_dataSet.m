% the dummy experiment has images files for each of the distinct conditions, 
% multiple wells, time points and color channels. 
% Each well is associated a specific file (e.g. listing a segmentation mask, etc)
nTime = 3;
nChannels = 2;
wells = {'A01','B08','D12'};
condition_names = {'control','drug 1','drug 2'};

%% generate folder
dirName = fullfile(pwd,'example_folder');
if ~exist(dirName)
    mkdir(dirName);
end
disp('********************************************************************');
disp(['Creating example folder holding dummy example files: ',dirName]);
curStr = 'Experiment has the following condition(s): ';
for i=1:numel(condition_names)
    curStr = [curStr,' ',condition_names{i}];
    if i<numel(condition_names)
        curStr = [curStr,' |'];
    end
end
disp(curStr);
curStr = 'the following well(s): ';
for i=1:numel(wells)
    curStr = [curStr,' ',wells{i}];
    if i<numel(wells)
        curStr = [curStr,' |'];
    end
end
disp(curStr);
disp(['Experiment has ',num2str(nTime),' time points and ',...
    num2str(nChannels),' channels.']);

%% generate dummy data files
for i=1:numel(condition_names)
    curDir = fullfile(dirName,condition_names{i});
    if ~exist(curDir,'dir')
        mkdir(curDir);
    end
    %generate ntime timpoints
    for j=1:nTime
        for k=1:numel(wells)
            for l=1:nChannels
                curStr = ['example file, condition ',condition_names{i},...
                    '; time ', num2str(j),'; well ',wells{k},'; channel ',num2str(l)];
                curFname = ['img',wells{k},num2str(l),'_t',num2str(j),'.txt'];
                fName = fullfile(curDir,curFname);
                fid = fopen(fName,'w');
                fprintf(fid, curStr);
                fclose(fid);
            end
            curStr = ['example mask, condition ',condition_names{i},...
                    '; well ',wells{k}];
            curFname = ['mask',wells{k},'.txt'];
            fName = fullfile(curDir,curFname);
            fid = fopen(fName,'w');
            fprintf(fid, curStr);
            fclose(fid);
        end 
    end
end

%% generate .ini config file

cfg = IniConfig();
cfg.AddSections({'Files','FileOptions','Output'});

% add the name pattern for Img file type
fpattern = [dirName,filesep,'{Condition}/img{Well #3}{Channel}_{Time}.txt'];
disp('"Image" files follow the pattern: ');
disp(['  ',fpattern]);
cfg.AddKeys('Files', {'Img'}, {fpattern});

% add the name pattern for Mask file type
fpattern = [dirName,filesep,'{Condition}/mask{Well #3}.txt'];
disp('"Mask" files follow the pattern: ');
disp(['  ',fpattern]);
cfg.AddKeys('Files', {'Mask'}, {fpattern});

% add the output folder
cfg.AddKeys('Output', {'outFolder'},{dirName});

% add the options
cfg.AddKeys('FileOptions', {'batchMode','recursive','verbose','useTerminal'}, ...
    {1,1,0,1});

% write to file
config_fname = fullfile(dirName,'cfg.ini');
cfg.WriteFile(config_fname);

disp('Created dummy example dataset and corresponding config file: ');
disp(['  ',fullfile(dirName,'cfg.ini')]);
disp(' ');
%% de-clutter workspace
clear i j k l curFname fid curStr fName nChannels nTime wells 
clear cfg condition_names fpattern curDir dirName