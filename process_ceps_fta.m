function varargout = process_ceps_fta(varargin )
% process_cps_fta: Computes the PSD of Frequency-Tagged data segments 
%                   marked by "Boundary" (or anyother user-defined label)
%                   in a moving-window fashion, and averages them.
%            
%
% @=============================================================================
% This function is a custom Brainstorm Process used to analyze
% Frequency-Tagged EEG data.
% 
%
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
%
% =============================================================================@
%
% Authors: Velu Prabhakar Kumaravel & Marco Buiatti
%
eval(macro_method);
end
%% ===== GET DESCRIPTION =====
function sProcess = GetDescription()
    % Description the process
    sProcess.Comment     = 'Compute Evoked Power Spectrum';
    sProcess.FileTag     = 'fta_cps';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = 'Frequency Tagging Analysis';
    sProcess.Index       = 600;
    sProcess.Description = 'placeholder';
    
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'results'};
    sProcess.OutputTypes = { 'timefreq','timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    
    % Options: Sensor types
    sProcess.options.sensortypes.Comment = 'Sensor Types or Names (empty=all): ';
    sProcess.options.sensortypes.Type    = 'text';
    sProcess.options.sensortypes.Value   = 'MEG, EEG';
    sProcess.options.sensortypes.InputTypes = {'data', 'results'};
    sProcess.options.sensortypes.Group   = 'input';
    

    % Options: Boundary event Label
    sProcess.options.BL.Comment = 'Boundary Event Label';
    sProcess.options.BL.Type    = 'text';
    sProcess.options.BL.Value   = 'boundary';

    % Options: Stimulation event Label
    sProcess.options.SL.Comment = 'Stimulation Event Label';
    sProcess.options.SL.Type    = 'text';
    sProcess.options.SL.Value   = 'DIN1'; % {Other options: 'DIN2', 'DIN3'}

    % === Time window
    sProcess.options.windowlen.Comment = 'Window Length';
    sProcess.options.windowlen.Type    = 'value';
    sProcess.options.windowlen.Value   = {10, 'seconds', 3};
    
    % === Upper bound overlap framing
    sProcess.options.Uepochoverlap.Comment = 'Max. Overlap Factor';
    sProcess.options.Uepochoverlap.Type    = 'value';
    sProcess.options.Uepochoverlap.Value   = {75, '%', 0};
    
    % === Lower bound overlap framing
    sProcess.options.Lepochoverlap.Comment = 'Min. Overlap Factor';
    sProcess.options.Lepochoverlap.Type    = 'value';
    sProcess.options.Lepochoverlap.Value   = {50, '%', 0};
    
    % === Zero-Padding Checkbox
    sProcess.options.isZeroPad.Comment     = 'Zero-Padding (in case of shorter epochs)';
    sProcess.options.isZeroPad.Type        = 'checkbox';
    sProcess.options.isZeroPad.Value       = 0;
    sProcess.options.isZeroPad.Controller  = 'cZeroPad';
    
    % === Zero-Padding
    sProcess.options.minlengthzeropad.Comment     = 'Minimum Epoch Length';
    sProcess.options.minlengthzeropad.Type        = 'value';
    sProcess.options.minlengthzeropad.Value       = {5, 'seconds', 1};
    sProcess.options.minlengthzeropad.Class       = 'cZeroPad';
    
end
%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) 
     
Comment = sProcess.Comment;

end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInput)
   % Initializing returned values
   OutputFiles = {};
   imagingKernel={};
   RDataMat={};
   Filepath=sInput.FileName;
   %Loading Data Matrix
   if (strcmp(sInput.FileType,'data'))
       DataStruct = in_bst_data(Filepath);
   elseif strcmp(sInput.FileType,'results')
       DataStruct = in_bst_data(sInput.DataFile);%original data struct
       %extracting result path
       [ ~ ,Filepath]=strtok(Filepath, '|');%removing "link" text from string
       [Rpath,Filepath]=strtok(Filepath, '|');%storing the right filepath of the result
       RDataMat=in_bst_data(Rpath);
       imagingKernel=RDataMat.ImagingKernel;
       Filepath=strtok(Filepath, '|');%path to the 'result' file, Rpath=path to the sensor file
   end
   
   inputData = DataStruct.F;
   
    % Compute Sampling Rate from the Time vector
    if (~isfield(DataStruct, 'Time') || isempty(DataStruct.Time))
        bst_report('Error', sProcess, sInput, 'No Time Points exist in the data');
        return;
    else
        dT  = mean(diff(DataStruct.Time));
        sRate  = round(1/dT);
        strMsg = ['The sampling rate is ' num2str(sRate) ''];
        bst_report('Info', sProcess, sInput, strMsg);
    end
    
    % Extracting data Boundaries
    if(~isempty(sProcess.options.BL.Value))
        user_event = sProcess.options.BL.Value;
    end

    % Extracting stimulation event label (cycle onset marker)
    stim_event = sProcess.options.SL.Value;
    if isempty(stim_event)
        bst_report('Error', sProcess, sInput, 'Stimulation Event Label cannot be empty.');
        return;
    end

    % Collect stimulation event sample indices
    [~, NEvents] = size(DataStruct.Events);
    stimSamples = [];
    isStim = 0;

    for StimPos = 1:NEvents
        if strcmp(stim_event, DataStruct.Events(StimPos).label)
            isStim = 1;

            if isfield(DataStruct.Events, 'samples') && ...
                    ~isempty(DataStruct.Events(StimPos).samples)

                stimSamples = DataStruct.Events(StimPos).samples;
            else
                stimSamples = round( ...
                    DataStruct.Events(StimPos).times * sRate);
            end

            % the first row
            stimSamples = sort(stimSamples(1, :));
            break
        end
    end

    if ~isStim
        strMsg = ['No ' stim_event ...
            ' event found in the data. Cannot align segments to the ' ...
            'stimulation cycle. Aborting...'];
        bst_report('Error', sProcess, sInput, strMsg);
        return;
    end

    isBound = 1;
    BoundPos = find(strcmp(user_event, {DataStruct.Events.label}), 1);
    if isempty(BoundPos)
        strMsg = ['No ' user_event ...
            ' event found in the data, processing all data points.'];
        bst_report('Warning', sProcess, [], strMsg);
        isBound = 0;
    end

    % Collect the raw (boundary-delimited) chunks
    if(isBound)
        if isfield(DataStruct.Events, 'samples') 
            seg_intervals=DataStruct.Events(BoundPos).samples;
        else
            seg_intervals=round(DataStruct.Events(BoundPos).times * sRate);
        end
        nSegments = size(seg_intervals, 2);
        rawChunks = cell(1, 0);

        for iSeg=1:nSegments-1
                start_idx = seg_intervals(1, iSeg);
                end_idx = seg_intervals(1, iSeg + 1) - 1;
                rawChunks{end+1}   = [start_idx, end_idx]; %#ok<AGROW>
        end
        
        if(seg_intervals(1) ~= 1)
             rawChunks{end+1} = [1, seg_intervals(1)];
        end
        
        if(seg_intervals(end) ~= size(inputData,2))
            rawChunks{end+1} = [seg_intervals(1, end), size(inputData,2)];
        end

        % Realign every chunk (inter-boundary AND the pre-first/post-last
        % edge pieces) to the first stimulation event it contains
        segDataStruct = cell(1, 0);
        for iChunk = 1:length(rawChunks)
            c_start = rawChunks{iChunk}(1);
            c_end   = rawChunks{iChunk}(2);
            firstStim = stimSamples(find(stimSamples >= c_start & stimSamples <= c_end, 1, 'first'));
            if isempty(firstStim)
                strMsg = sprintf(['No ' stim_event ' event found within chunk [%d %d] samples. ' ...
                    'Discarding this chunk.'], c_start, c_end);
                bst_report('Warning', sProcess, [], strMsg);
                continue;
            end
            segDataStruct{end+1} = inputData(:, firstStim:c_end); %#ok<AGROW>
        end
        
        inputData = segDataStruct; % override inputData in case of boundary segments
    else
        % No boundary events: still need to align the single chunk
        % (whole recording) to the first stimulation event
        c_start = 1;
        c_end = size(inputData, 2);
        firstStim = stimSamples(find(stimSamples >= c_start & stimSamples <= c_end, 1, 'first'));
        if isempty(firstStim)
            strMsg = ['No ' stim_event ' event found in the data. Cannot align segments. Aborting...'];
            bst_report('Error', sProcess, sInput, strMsg);
            return;
        end
        inputData = {inputData(:, firstStim:c_end)}; % else convert the data into a cell array
    end

    if(isempty(inputData))
       strMsg = 'No segment contained a stimulation event; nothing to process.';
       bst_report('Error', sProcess, [], strMsg);
       return;
    end
    
   %Window is positive and not =0
   WindowLength=sProcess.options.windowlen.Value{1}*sRate; % in samples
   if(WindowLength < 0)
       bst_report('Error', sProcess, sInput, 'Window Length cannot be zero or negative. Aborting...');
       return;
   end
    
   % Zero-padding
   if isfield(sProcess.options, 'isZeroPad') && isfield(sProcess.options.isZeroPad, 'Value') && ~isempty(sProcess.options.isZeroPad.Value)
       isZeroPad = sProcess.options.isZeroPad.Value;
   else
       isZeroPad = 1;
   end
   
   if(isZeroPad)
       min_length_zp = sProcess.options.minlengthzeropad.Value;
   else
       min_length_zp = 0;
   end
   
   inputData_copy = inputData; % use the copy to perform zero-padding
   
   
   if (isZeroPad) 
       inputData_pad = [];
       count = 1; % Count of zero-padded (new) segments
       for iSeg=1:length(inputData_copy)
           [ChannelNumber, Datasize]= size(inputData_copy{iSeg});
           if(Datasize < WindowLength && Datasize >= min_length_zp{1}*sRate && Datasize <= WindowLength)
               inputData_pad{count} = [inputData_copy{iSeg} zeros(ChannelNumber, WindowLength - size(inputData_copy{iSeg},2))]; %Zero-padding
               strMsg = ['Data Segment : ' num2str(iSeg) ' is zero-padded'];
               bst_report('Warning', sProcess, [], strMsg);
           else
               inputData_pad{count} = inputData_copy{iSeg};
           end
           count = count + 1;
       end
       inputData = inputData_pad;
   end
   
   % Logic to check if none of the segments meet the window length criteria
   for iSeg = 1:length(inputData)
       size_factor(iSeg)=length(inputData{iSeg}(1,:))/WindowLength;
   end
   
   inputData(size_factor<1) = [];
   
   if(isempty(inputData))
       strMsg = 'Too short data for the PSD computation';
       bst_report('Error', sProcess, [], strMsg);
       return;
   end
   
   % Taper is Square
   psdKernel   = ones(1,WindowLength)/WindowLength; 
   
      %lower limit<upper limit
   if (sProcess.options.Lepochoverlap.Value{1}>sProcess.options.Uepochoverlap.Value{1})
       bst_report('Warning', sProcess, sProcess.options.Lepochoverlap.Value{1}, 'Upper overlap bound lower than lower overlap bound, automatic switch');
        temp=sProcess.options.Lepochoverlap.Value{1};
        sProcess.options.Lepochoverlap.Value{1}=sProcess.options.Uepochoverlap.Value{1};
        sProcess.options.Uepochoverlap.Value{1}=temp;
   end
   
  %make the code easier to read
  UpperBound=sProcess.options.Uepochoverlap.Value{1}/100;
  LowerBound=sProcess.options.Lepochoverlap.Value{1}/100;

  % find bad channels (if any)
  badChannels = find(DataStruct.ChannelFlag==-1);
   
  [ps, interval, Nwin]=Compute(sProcess,sInput,inputData, WindowLength, ...
      UpperBound, LowerBound, psdKernel, imagingKernel, badChannels);
  
   %---Power Spectrum---
   ps = permute(ps, [1 3 2]); 
   %Output file creation
   f=interval*sRate/WindowLength;
   [ChannelNumber,~ ,~]=size(ps);
   Rows=1:1:ChannelNumber;
   FileMat = db_template('timefreqmat');
   OutputFiles=bst_process('GetNewFilename',bst_fileparts(Filepath),'timefreq_psd');
   FileMat.ChannelFlag=DataStruct.ChannelFlag;
   FileMat.TF=ps;
   FileMat.Comment= sprintf('PSD: %d/%d ms %s',sum(Nwin),WindowLength*(1/sRate)*1000,'Power');
   FileMat.DataType='data';
   FileMat.Time=[0,size(DataStruct.F,2)*(1/sRate)];
   FileMat.Freqs=f;
   FileMat.RowNames=Rows;

   if(strcmp(sInput.FileType,'data')||strcmp(sInput.FileType,'raw')) % sensor-level analysis
       ChannelFile = bst_get('ChannelFileForStudy', sInput.FileName);
       if isempty(ChannelFile)
           error('No channel definition available for this file.');
       end
       % Load channel file
       ChannelMat  = in_bst_channel(ChannelFile);
       % Get channels we want to process
       iChannels = 1:length(ChannelMat.Channel);
       FileMat.RowNames = {ChannelMat.Channel(iChannels).Name};
   end
   
   FileMat.Measure='power';
   FileMat.Method='psd';
   FileMat.DataFile=sInput.FileName;
   FileMat.nAvg=DataStruct.nAvg;
   FileMat.Options=sProcess.options;
   FileMat.History=DataStruct.History;
   if(~isempty(RDataMat))
       FileMat.DataType='results';
       FileMat.ChannelFlag=RDataMat.ChannelFlag;
       FileMat.Function=RDataMat.Function;
       FileMat.GoodChannel=RDataMat.GoodChannel;
       FileMat.SurfaceFile=RDataMat.SurfaceFile;
       FileMat.HeadModelFile=RDataMat.HeadModelFile;
       FileMat.HeadModelType=RDataMat.HeadModelType;
       FileMat.Whitener=RDataMat.Whitener;
       FileMat.History=RDataMat.History;
   end
   FileMat = bst_history('add', FileMat, 'compute', 'Frequency Tagging Analysis - Power Spectrum');
   % Save the new file
    save(OutputFiles,'-struct','FileMat');
    % Reference OutputFile in the database:
    db_add_data(sInput.iStudy, OutputFiles,FileMat);
end

%% ===== COMPUTE =====
function [ps, interval, Nwin]=Compute(sProcess,sInput,inputData, ...
    WindowLength, UpperBound, LowerBound, psdKernel, imagingKernel, badCh)
    
    [ps, interval, Nwin] = fta_ps(inputData, WindowLength, UpperBound, ...
        LowerBound, psdKernel, imagingKernel, badCh);    
    bst_report('Info', sProcess, sInput, sprintf('Number of windows used: %d\n',Nwin));
    
end


%% ===== PSD Adaptive Windowing Logic goes here =====

function [ps, interval, nWin] = fta_ps(inputData, winLen, uppBound, lowBound, psdKernel, imagingKernel, badCh)
 
    % imagingKernel(:, badCh) = [];
    for iSeg = 1:length(inputData)
        inputData{iSeg}(badCh, :) = []; 
    end
    
    [ChNumber, ~]= size(inputData{1});
    nWin = zeros(1,length(inputData));

    w = winLen*(sum(psdKernel.^2));	%window squared and summed
    %% EPS %%
    tr = 0;
    for iSeg = 1:length(inputData)
        DataLength = size(inputData{iSeg}(1,:),2);
        nWin(iSeg)=floor(DataLength/winLen);
        for i=1:nWin(iSeg)
            tr=tr+1;
            for ch = 1:ChNumber
                wd(ch,:,tr)=psdKernel.*inputData{iSeg}(ch,(i-1)*winLen+1:(i*winLen));
            end
        end
        bst_progress('inc',  ceil(1/length(inputData)*100));
    end
    wd_mean=mean(wd,3);
    fwd_mean=fft(wd_mean,[],2);
    if ~isempty(imagingKernel)
        fwd_mean = imagingKernel * fwd_mean; % computing the source psd
    end
    ap=abs(fwd_mean).^2;
    n=1;
    
    if mod(winLen, 2) == 0
        ps(:,1) = ap(:,1) / (w * n);
        ps(:, 2:winLen/2) = 2*ap(:, 2:winLen/2) / (w * n);
        ps(:, (winLen/2) + 1) = ap(:, (winLen/2) + 1) / (w * n);
        interval = 0:1:(winLen/2);
    else
        ps(:,1) = ap(:,1) / (w * n);
        ps(:, 2:ceil(winLen/2)) = 2*ap(:, 2:ceil(winLen/2)) / (w * n);
        ps(:, ceil(winLen/2) + 1) = ap(:, ceil(winLen/2) + 1) / (w * n);
        interval = 0:1:(ceil(winLen/2));
    end
end