%% Armeo_2DHorizCatch
%
% This script allows you to manage the raw data from .csv files generated 
% by the 'Armeocontrol 1.22' software provided with the ArmeoSpring.
% It was designed to import and analyze individual exercises (trials)
% composing an evaluation session with the '2D-horizontal catch' exercise.
%
% The "data" dataset offered with the script allows to browse csv file 
% content and to test the script
%
%%%% Details concerning the '2D-horizontal catch' exercise %%%%%%%%%%%%%%%%
% - The 2D-horizontal pointing task required to move the cursor in order to 
% catch balls that appeared sequentially on the screen.
% - When a ball was caught, it disappeared and another appeared at a new 
% fixed location
% - During a trial, 12 balls had to be caught and the 
% time to catch a ball was limited to 10 seconds. If it was exceeded, the 
% ball disappeared and another ball appeared at the new location.
% - Only the last 11 targets (balls) are considered for analysis.
% - Trajectory to catch a ball is computed only when the preceding ball was 
% successfully caught
%
% ****************************************************************************
% To launch the procedure, you just have to run the script 'Armeo_2DHorizCatch'
% *****************************************************************************
%
% A synthesis table containing the computed parameters for each trial (one
% csv file per trial) is generated in txt format, and saved in the 'export'
% sub-folder, itself contained in the same folder as the script
% If the 'export' folder does not exist, it will be generated automatically 
%
% During the script process, you will be asked to whether you want to:
%   - Graphically visualize 2D movements' trajectories used for the 
%   calculation of kinematic parameters
%   - Save plots of these 2D trajectories (in the 'export' folder)
%
%   This script has been developed within the framework of the publication 
%   of the article: N. Brihmat, I. Loubinoux, E. Castel-Lacanal, P. Marque,
%   D. Gasq. Kinematic parameters obtained with the ArmeoSpring for 
%   upper-limb assessment after stroke: A reliability and learning effect 
%   study guiding the parameters' use. J. NeuroEng Rehab. 2020
%
%   The objective is to facilitate the use and analysis of raw ArmeoSpring 
%   data in clinical practice as proposed in the paper.
%
% Authors: David Gasq[1,2], MD, PhD & Nabila Brihmat[1], PhD - May 2020
% [1] ToNIC, Toulouse NeuroImaging Center, Université de Toulouse, Inserm, UPS, France
% [2] University Hospital of Toulouse, Hôpital de Rangueil, Toulouse, France
% email contact: gasq.d@chu-toulouse.fr
% 
% Copyright (C) 2020  David Gasq
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

%% -1- Path management
[current_path, ~, ~] = fileparts(which('Armeo_2DHorizCatch.m'));
settings.current_path = current_path;
cd(current_path)
exportdir = 'export';
if exist(exportdir,'dir') == 0
    cd(current_path)
    mkdir export
end
settings.export_path = sprintf('%s%s%s',current_path,filesep,exportdir);

%% -2- Run Armeo_kinematics
[Armeo_Kin] = Armeo_kinematics(settings);

clear exportdir current_path settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Armeo_Kin] = Armeo_kinematics(settings)

% Function Armeo_kinematics
% This function import raw data from .csv files generated 
% by the 'Armeocontrol 1.22' software provided with the ArmeoSpring.
% It was designed to import and analyze individual exercises (trials)
% composing an evaluation session with the '2D-horizontal catch' exercise.
%
% Inputs:  
% - 'settings' is a structure with the following fields: 
%   - current_path: path to this script
%   - export_path: path to the folder where the csv file containing the computed
%   data will be saved
%
% Output:
%   - Armeo_Kin: a table with the computed variables for each trial. 
%       Variables are the mean of successfully caught balls accross the task: 
%       number of success, TotalTaskTime (s), Score (%), MovementTime (s), 
%       PeakVel (cm/s), HandPathRatio, nVelPeaks
%   The table is saved in .txt in the folder 'export'
%   - Inspection plots of 2D trajectories from which the parameters have 
%       been computed for each trial (1 plot per csv file processed)are
%       available for simple visual inspection and/or saving in .png format
%
%   Users can thus use the variables they wish, by averaging the variables 
%   on the trials of their choice.
%
%   Author: David Gasq, MD, PhD
%   ToNIC, Toulouse NeuroImaging Center, Université de Toulouse, Inserm, UPS, France
%   University Hospital of Toulouse, Hôpital de Rangueil, Toulouse, France
%   email address: gasq.d@chu-toulouse.fr
%   October 2018; Last revision: May 2020

%% Path & settings
current_path = settings.current_path ;
export_path = settings.export_path ;
titlebox = 'Select one or more csv files';
[files,import_path,index] = uigetfile('.csv','title',titlebox,...
    'MultiSelect','on') ;
if isequal(files,0) || isequal(index,0) || isequal(files,titlebox)
   errordlg('No csv file has been selected: Please start again',...
       'File selection error');
   Armeo_Kin = NaN;
   return
else
    if ~ iscell(files)
        files = {files};
    end
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    fprintf('%d csv files will be processed \n',numel(files));
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp(' ');
end
choice = questdlg(sprintf('%s %s',...
    'Would you like to plot and save 2D trajectories',...
    'processed for each trial (one trial = 1 csv file)?'),...
'Inspection plot','Just plot','Plot and save','No','Plot and save');
switch choice
    case 'Just plot'
        inspectPlot = 1;
    case 'Plot and save'
        inspectPlot = 2;
    case 'No'
        inspectPlot = 0;
end

%% Data processing
data_files = cell(1,numel(files));
for n = 1:1:numel(files)

    filename = sprintf('%s',files{1,n});
    fprintf('"%s" is being processed \n',filename);
    filenamepath = fullfile(import_path,filename);
    
    data_files{1,n} = armeo2D_process(filenamepath);

    clear filename filenamepath

end % for n = 1:1:numel(files)

%% Graphical inpection
if inspectPlot == 1 || inspectPlot == 2
    
    % settings
    %list_color_name = {'Blue', 'Magenta', 'DarkGreen','Grey','Cyan',...
    %               'Orange','Black','Green','Purple','Brown','Violet'} ;
    list_color_rgb = {[0 0 1],[1 0 1],[0 0.39 0],[0.5 0.5 0.5],[0 1 1],...
                [1 0.65 0],[0 0 0],[0 0.5 0],[0.5 0 0.5],[0.65 0.16 0.16],...
                [0.93 0.510 0.93]} ;
    list_objects = {'Object 2','Object 3','Object 4','Object 5',...
                     'Object 6','Object 7','Object 8','Object 9',...
                     'Object 10','Object 11','Object 12'} ;

    sp = 11; % size of police
    fn = 'Times-Roman';
    
    % Loop over trials
    cd(export_path)
    for nt = 1:1:numel(data_files)

        plotArmeo2D = figure('NumberTitle','off',...
            'Name',sprintf('HorizCatch_trial%d_%s',nt,files{1,n}),...
            'units','normalized',...
            'outerposition',[0 0 1 1]);
        hold on;

        list_label = {};
        % loop over objects
        for nob = 1:1:numel(data_files{1, nt}.XZ)

            if ~isnan(data_files{1, nt}.XZ{1,nob})
                X = data_files{1, nt}.XZ{1,nob}(:,1);
                Z = data_files{1, nt}.XZ{1,nob}(:,2);
                plot(X,Z,'color',list_color_rgb{nob});
                list_label{numel(list_label)+1} = list_objects{1,nob} ;
                clear X Z
            end % if ~isnan(D.XY{ne,nob} == 1   

        end % for nt = 1:1:numel(data_files{1, 1}.XZ)

        xlabel('Medio-lateral axis (cm)','Color','k','FontSize',sp,'FontWeight','bold','FontName',fn);
        ylabel('Proximo-distal axis (cm)','Color','k','FontSize',sp,'FontWeight','bold','FontName',fn);  
        title(sprintf('Horizontal Catch trial %d (%s)',nt,files{1,n}),...
            'Interpreter', 'none');
        legend(list_label,'Location','best','Orientation','vertical');
        
        if inspectPlot == 2
            NamePlot = sprintf('Plot_HorizCatch_%s_trial%d',files{1,n},nt);
            NamePlot = strrep(NamePlot,'.csv','');
            print(plotArmeo2D,NamePlot,'-dpng','-r300');
        end % if inspectPlot == 2
        
        uicontrol('Parent',plotArmeo2D,'Style','pushbutton',...
            'String','Click here to close the figure and continue',...
            'BackgroundColor','green','ForegroundColor','black',...
            'FontWeight','bold','FontSize',14,'Units','normalized',...
            'Position',[0.6 0.95 0.3 0.05],'callback','close(''gcbf'')');
        waitfor(plotArmeo2D);
        close all
        clear list_label
                
    end  % for nt = 1:1:numel(data_files)
    cd(current_path)

end % if inspectPlot == 1 || inspectPlot == 2

%% Synthesis of results
varNames = {'numTrial' 'n_objects' 'n_success' ...
                'TotalTaskTime' 'Score' 'MovementTime' 'PeakVel' ...
                'HandPathRatio' 'nVelPeaks'} ;

% loop over trials
data = nan(numel(data_files),numel(varNames)) ;
for nt = 1:1:numel(data_files)
    
    data(nt,1) = nt;
    data(nt,2) = 11;
    data(nt,3) = nansum(data_files{1, nt}.failSucc);
    data(nt,4) = data_files{1, nt}.totalDuration ;
    data(nt,5) = data_files{1, nt}.score ;
    data(nt,6:9) = nanmean(data_files{1,nt}.data(:,2:5),1) ;
    
end % for nt = 1:1:numel(data_files)

% Table of computed data
Armeo_Kin = array2table(data) ;
Armeo_Kin.Properties.VariableNames = varNames ;
Armeo_Kin.Properties.RowNames = files;

cd(export_path)
d = datetime('now','Format','yyyMMddHHmmss');
writetable(Armeo_Kin,sprintf('Armeo2D_%s.txt',d),...
            'Delimiter',',',...
            'WriteRowNames',true);

end % function [Armeo_Kin] = Armeo_kinematics(settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function D = armeo2D_process(filename)
    
%   Function armeo2D_process
%   Read and process Armeo data from a csv file
%
%   Inputs:
%   - filename: 'path/filename' of the csv file (str)
%
%   Output:
%   - D: structure with processed data
%       - D.no: targets' number
%       - D.totalDuration: total duration of exercise (11 last targets)
%       - D.ObjDuration: duration for each object caught (11 last targets),
%           even if the catch of the preceding object fails 
%       - D.score: Armeo's score considering the 11 last targets
%       - D.failSucc: failure (0) or success (1) for each objet
%       - D.label: label of computed kinematic parameters
%       - D.data: matrix of computed kinematic parameters for each object,
%           computed only when the preceding object was successful caught
%
%   Author: David Gasq, MD, PhD
%   ToNIC, Toulouse NeuroImaging Center, Université de Toulouse, Inserm, UPS, France
%   University Hospital of Toulouse, Hôpital de Rangueil, Toulouse, France
%   email address: gasq.d@chu-toulouse.fr
%   October 2018; Last revision: May 2020

    %% -1- read data from csv file
    Arm_2Ddat = readArmeo2Ddata(filename);
    TS = Arm_2Ddat.Timestamp;
    time = nan(numel(TS),4);

    for i = 1:1:numel(TS) 
        TSi = TS{i};
        H = str2double(TSi(1:2));
        M = str2double(TSi(4:5));
        S = str2double(TSi(7:12));
        tictac = H*3600 + M*60 + S;
        time(i,:) = [H M S tictac];
        clear TSi H M S tictac
    end

    Difftime = [diff(time(:,4)); NaN] ;
    fs = round(1/nanmean(Difftime),0);

    % Raw data xz
    XZ = [Arm_2Ddat.x Arm_2Ddat.z];

    % Index object apparition
    ObjApp = Arm_2Ddat.Start;
    ObjAppIdx = find(ObjApp == 1);
    IdxApp = [1;ObjAppIdx;NaN];

    % Index Objectcaught & Objectdisappeared
    ObjCau = Arm_2Ddat.Objectcaught;
    ObjCauIdx = find(ObjCau == 1);

    ObjDis = Arm_2Ddat.Objectdisappeared;
    ObjDisIdx = find(ObjDis == 1);

    ObjCD = sort([ObjCauIdx;ObjDisIdx],1,'ascend');
    IdxCauDis = [NaN;ObjCD];

    CauDis10 = nan(numel(IdxCauDis),1);
    for i = 1:1:numel(IdxCauDis) 
        if ismember(IdxCauDis(i),ObjCauIdx)
            CauDis10(i,1) = 1;
        elseif ismember(IdxCauDis(i),ObjDisIdx)
            CauDis10(i,1) = 0;
        else
            CauDis10(i,1) = NaN;
        end  
    end

    not = numel(IdxApp(~isnan(IdxApp))); % 12 targets
    if not == 12
        no = not - 1; % The 11 last targets are considered for analysis
    else
        sprintf('Error in PP_Arm3D because because there is not 12 targets')
        return
    end

    %% -2- data processing
    D.no = no ;
    % extraction of data corresponding to the last 11 targets
    totalDuration = size(XZ(IdxApp(2):IdxCauDis(13),:),1)/fs; 
    D.totalDuration = totalDuration;
    ObjDuration = nan(1,no) ;

    success = nansum(CauDis10(3:13)) ;
    failure = no - success ;
    score = ((no - failure) / no )*100 ;
    D.score = score ;
    failSucc = CauDis10(3:13) ;
    D.failSucc = failSucc ;
    
    for nn = 1:1:no

        Ini = IdxApp(nn+1);
        End = IdxCauDis(nn+2);
        CD10 = CauDis10(nn+2);
        CD10past = CauDis10(nn+1);
        failSucc_n = failSucc(nn) ;
        vXZ = XZ(Ini:End,:);
        ObjDuration(1,nn) = size(vXZ,1)/fs ;

        if nn == 1 % first iteration

            if CD10 == 1 % success

                P = hand2D(vXZ,fs,failSucc_n);
                if nn == 1
                    D.label = P.label ;
                end
                D.data(nn,1:numel(P.data)) = P.data ;
                D.XZ{1,nn} = vXZ ;

            elseif CD10 == 0 % failure

                P = hand2D(nan,fs,failSucc_n);
                if nn == 1
                    D.label = P.label ;
                end
                D.data(nn,1:numel(P.data)) = P.data ;
                D.XZ{1,nn} = nan ;

            else
                error('Error in CD10 value for nn==1 in "armeo2D_process"');
            end % if CD10 == 1

        elseif nn > 1

           if CD10 == 1 && CD10past == 1 % success

                P = hand2D(vXZ,fs,failSucc_n);
                D.data(nn,1:numel(P.data)) = P.data ;
                D.XZ{1,nn} = vXZ ;

            elseif CD10 == 0 || CD10past == 0 % failure

                P = hand2D(nan,fs,failSucc_n);
                D.data(nn,1:numel(P.data)) = P.data ;
                D.XZ{1,nn} = nan ;

            else
                error('Error in CD10 value for nn>1 in "armeo2D_process"');
            end % if CD10 == 1 

        end % if nn == 1
        clear Ini End CD10 vXZ

    end % for nn = 1:1:no

    D.ObjDuration = ObjDuration ;

end % D = armeo2D_process(filename, ...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T = readArmeo2Ddata(filename)
%   Function readArmeo2Ddata
%   Import numeric data from a text file as a matrix
%   "T = readArmeo2Ddata(filename)": Reads data from text file 'filename' 
%   
%   Inputs:
%   - filename: 'path/filename' of the csv file (str)
%
%   Output:
%   - T: table with data from the csv file

%   Examples:
%   T = readArmeo2Ddata('path.../AE1_Hor_20180525132838_Data.csv');
% 
%   Author: David Gasq, MD, PhD
%   ToNIC, Toulouse NeuroImaging Center, Université de Toulouse, Inserm, UPS, France
%   University Hospital of Toulouse, Hôpital de Rangueil, Toulouse, France
%   email address: gasq.d@chu-toulouse.fr
%   May 2019; Last revision: June 2020

    try
        delimiterIn = ';';
        headerlinesIn = 1;
        M = importdata(filename,delimiterIn,headerlinesIn);
        T = table;
        T.Timestamp = M.textdata(2:end);
        T.x = M.data(:,1);
        T.z = M.data(:,2);
        T.Objectnumber = M.data(:,3);
        T.Start = M.data(:,4);
        T.Objectcaught = M.data(:,5);
        T.Objectdisappeared = M.data(:,6);
    
    catch
        
        error('Error in reading the csv file');
        T = NaN;
        return
    
    end % try
    
end % function T = readArmeo2Ddata(...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function P = hand2D(vXZ,fs,failSucc_n)

%   Function hand2D
%   Computation of kinematic parameters from 2D Armeo raw data
%     
%   Inputs:
%   - vXZ = 2D vector if success, else Nan
%   - fs: frequency of sampling
%   - failSucc_n = failure (0) or success (1 = catch)
%
%   Output:
%   - P.label: label of computed kinematic parameters
%   - P.data: matrix of computed kinematic parameters
%
%   Author: David Gasq, MD, PhD
%   ToNIC, Toulouse NeuroImaging Center, Université de Toulouse, Inserm, UPS, France
%   University Hospital of Toulouse, Hôpital de Rangueil, Toulouse, France
%   email address: gasq.d@chu-toulouse.fr
%   May 2019; Last revision: May 2020

if ~isnan(vXZ) % success
    
    duration = size(vXZ,1)/fs;
    
    % Lowpass filter: order = 2; cut frequency = 6Hz
    [b,a]=butter(2,(6*2)/fs,'low');
    vXZf = filtfilt(b,a,vXZ);
        
    % 2D scalar length
    Lmin = sqrt((vXZf(end,1)-vXZf(1,1))^2+(vXZf(end,2)-vXZf(1,2))^2);
    % 2D covered length
    Lcov = sum(sqrt(diff(vXZf(:,1)).^2+diff(vXZf(:,2)).^2));

    % PathRatio
    PathRatio = Lcov/Lmin ; % expressed from 1 to infinity

    % 2D scalar velocity, nVelPeak and PeakVel
    velf = nan(size(vXZf,1),1);
    for p = 1:1:size(vXZf,1)-1
        velf(p,1) = sqrt((vXZf(p+1,1)-vXZf(p,1))^2+(vXZf(p+1,2)-vXZf(p,2))^2)/(1/fs);
    end
    
    [nVelPeaks,~,maxVel] = findPpeaks(velf);

    P.label = {'SuccessFailure',...
            'MovementTime',...
            'PeakVel',...
            'HandPathRatio',...
            'nVelPeaks'};
    P.data = [failSucc_n,...
            duration,...
            maxVel,...
            PathRatio,...
            nVelPeaks];

else % failure
    
    P.label = {'SuccessFailure',...
            'MovementTime',...
            'PeakVel',...
            'HandPathRatio',...
            'nVelPeaks'};
    P.data = [failSucc_n,...
            nan,...
            nan,...
            nan,...
            nan];
    
end

end % function P = hand2D(vXZ,...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [nPPks,PPksVal,maxPPksVal] = findPpeaks(V)

%   Function findPpeaks
%   Extraction of positive peaks from a signal
%     
%   Inputs:
%   - V = one-row data vector
%
%   Output:
%   - nPPks: numbers of positive peaks
%   - PPksVal: positive peaks' values
%   - maxPPksVal: maximal positive peak's value
%
%   Author: David Gasq, MD, PhD
%   ToNIC, Toulouse NeuroImaging Center, Université de Toulouse, Inserm, UPS, France
%   University Hospital of Toulouse, Hôpital de Rangueil, Toulouse, France
%   email address: gasq.d@chu-toulouse.fr
%   July 2013; Last revision: May 2020

if size(V,1) > 1
    Vp = V';
else
    Vp = V;
end

sa = sign(diff([-inf Vp]));
PP = find(diff(sa) == -2);

% Remove last and first points of the curve if considered as a peak 
if ~isempty(PP)
    if PP(end) == length(Vp)
        LocsPP = PP(1:end-1);
    else
        LocsPP = PP;
    end
    if PP(1) == 1
        LocsPP = LocsPP(2:end);
    end  
end

nPPks = length(LocsPP);
PPksVal = Vp(LocsPP);
maxPPksVal = max(PPksVal);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%