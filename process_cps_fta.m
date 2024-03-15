function varargout = process_cps_fta(varargin )
% process_cps_fta: Computes the PSD of Frequency-Tagged data segments marked by "Boundary" (or anyother user-defined label)
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
    sProcess.Comment     = 'Compute Power Spectrum';
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
    sProcess.options.BL.Comment = 'Event Label';
    sProcess.options.BL.Type    = 'text';
    sProcess.options.BL.Value   = 'boundary';
    
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
    
    %---Extracting data Boundaries
    
    user_event = sProcess.options.BL.Value;
    
    % Find the index of user_event in DataStruct.Events
    BoundPos = find(strcmp(user_event, {DataStruct.Events.label}), 1);
    
    % Check if user_event was not found
    if isempty(BoundPos)
        strMsg = ['No ''' user_event ''' event found in the data, processing all data points.'];
        bst_report('Warning', sProcess, [], strMsg);
        isBound = 0;
    else
        isBound = 1;
    end

    
    % Segment the data
    if(isBound)
        if isfield(DataStruct.Events, 'samples') % check with Marco
            seg_intervals=DataStruct.Events(BoundPos).samples;
        else
            seg_intervals=round(DataStruct.Events(BoundPos).times * sRate);
        end
        [~,nSegments]=size(seg_intervals);

        for iSeg=1:nSegments-1
                segDataStruct{iSeg}=inputData(1:size(inputData,1), (seg_intervals(iSeg)+1):(seg_intervals(iSeg+1)));
        end
        
        if(seg_intervals(1) ~= 1)
             segDataStruct{length(segDataStruct)+1}=inputData(1:size(inputData,1), 1:seg_intervals(1));
        end
        
        if(seg_intervals(end) ~= size(inputData,2))
            segDataStruct{length(segDataStruct)+1}=inputData(1:size(inputData,1), seg_intervals(end):size(inputData,2));
        end
        
        inputData = segDataStruct; % override inputData in case of boundary segments
    else
        inputData = {inputData}; % else convert the data into a cell array
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
   
  [ps, interval, Nwin]=Compute(sProcess,sInput,inputData, WindowLength, UpperBound, LowerBound, psdKernel, imagingKernel);
  
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
function [ps, interval, Nwin]=Compute(sProcess,sInput,inputData, WindowLength, UpperBound, LowerBound, psdKernel, imagingKernel)
    
    [ps, interval, Nwin] = fta_ps(inputData, WindowLength, UpperBound, LowerBound, psdKernel, imagingKernel);    
    bst_report('Info', sProcess, sInput, sprintf('Number of windows used: %d\n',Nwin));
    
end


%% ===== PSD Adaptive Windowing Logic goes here =====

function [ps, interval, nWin] = fta_ps(inputData, winLen, uppBound, lowBound, psdKernel, imagingKernel)

    [ChNumber, ~]= size(inputData{1});
    nWin = zeros(1,length(inputData));

    w = winLen*(sum(psdKernel.^2));	%window squared and summed

    ap = 0;
    for iSeg = 1:length(inputData)
        DataLength = size(inputData{iSeg}(1,:),2);
        if DataLength-winLen < floor(winLen*(1-uppBound))
            Nmax = 1; % Set to 1 for non-overlapping windows
            winStart = 0;
        else
            if uppBound == 0 && lowBound == 0
                % Set default values when both uppBound and lowBound are zero
                Nmax = 1; % Set to 1 for non-overlapping windows
                winStart = 0;
            else
                Nmax = max(1, ceil(DataLength / (winLen * max(eps, lowBound))));
                winStart = floor((DataLength - winLen) / max(1, (Nmax - 2)));
            end
        end


        nWin(iSeg)=max(1, nWin(iSeg)); % number of consecutive full windows

        ap_loc=0; %average power
        for i=1:nWin(iSeg)
            for ch = 1:ChNumber
                wd(ch,:)=psdKernel.*inputData{iSeg}(ch,(i-1)*winStart+1:(i-1)*winStart+winLen);
            end
            fwd = fft(wd,[],2);
            if ~isempty(imagingKernel)
                fwd = imagingKernel * fwd; % computing the source psd
            end
            ap_loc = ap_loc + abs(fwd).^2;
        end
        ap = ap + ap_loc; % output average power
        bst_progress('inc',  ceil(1/length(inputData)*100));
    end

    n=sum(nWin);
    
    if mod(winLen, 2) == 0
        ps(:,1) = ap(:,1) / (w * n);
        ps(:, 2:winLen/2) = (ap(:, 2:winLen/2) + ap(:, winLen:-1:(winLen/2) + 2)) / (w * n);
        ps(:, (winLen/2) + 1) = ap(:, (winLen/2) + 1) / (w * n);
        interval = 0:1:(winLen/2);
    else
        ps(:,1) = ap(:,1) / (w * n);
        ps(:, 2:ceil(winLen/2)) = (ap(:, 2:ceil(winLen/2)) + ap(:, winLen:-1:ceil(winLen/2) + 1)) / (w * n);
        ps(:, ceil(winLen/2) + 1) = ap(:, ceil(winLen/2) + 1) / (w * n);
        interval = 0:1:(ceil(winLen/2));
    end
end

function [ps, interval, nWin] = fta_ps_old(inputData, winLen, uppBound, lowBound, psdKernel, imagingKernel)

    [ChNumber, ~]= size(inputData{1});
    nWin = zeros(1,length(inputData));

    w = winLen*(sum(psdKernel.^2));	%window squared and summed

    ap = 0;
    for iSeg = 1:length(inputData)
        DataLength = size(inputData{iSeg}(1,:),2);
        if DataLength-winLen<floor(winLen*(1-uppBound))  % No overlapping if the data length is lesser than
            DataLength=winLen;
            Nmax=ceil(DataLength/(winLen*lowBound)); % max number of consecutive HALF windows (non-overlapping if l is multiple of wl/2)
            winStart=0;
        else
            Nmax=ceil(DataLength/(winLen*lowBound)); % max number of consecutive HALF windows (non-overlapping if l is multiple of wl/2)
            winStart=floor((DataLength-winLen)/(Nmax-2)); % step
        end

        nWin(iSeg)=Nmax-1; % number of consecutive full windows

        ap_loc=0; %average power
        for i=1:nWin(iSeg)
            for ch = 1:ChNumber
                wd(ch,:)=psdKernel.*inputData{iSeg}(ch,(i-1)*winStart+1:(i-1)*winStart+winLen);
            end
            fwd = fft(wd,[],2);
            if ~isempty(imagingKernel)
                fwd = imagingKernel * fwd; % computing the source psd
            end
            ap_loc = ap_loc + abs(fwd).^2;
        end
        ap = ap + ap_loc; % output average power
        bst_progress('inc',  ceil(1/length(inputData)*100));
    end

    n=sum(nWin);
    
    ps(:,1)=ap(:,1)/(w*n);
    ps(:,2:winLen/2)=(ap(:,2:winLen/2)+ap(:,winLen:-1:winLen/2 +2))/(w*n);
    ps(:,winLen/2 +1)=ap(:,winLen/2 +1)/(w*n);
    interval=0:1:(winLen/2);
end