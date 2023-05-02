function varargout = process_citc_fta(varargin )
% process_citc_fta: Computes Inter-Trial Coherence on all channels of EEG data over
% consecutive non-overlapping windows (trials)
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
    sProcess.Comment     = 'Compute Inter-Trial Phase Coherence';
    sProcess.FileTag     = 'citc';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = 'Frequency Tagging Analysis';
    sProcess.Index       = 603;
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

    % Options: Window Length
    sProcess.options.windowlen.Comment = 'Window Length';
    sProcess.options.windowlen.Type    = 'value';
    sProcess.options.windowlen.Value   = {5, 'seconds', 3};

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
    if(~isempty(sProcess.options.BL.Value))
        user_event = sProcess.options.BL.Value;
    end
    
    % Please optimize this code - this looping isn't necessary!
    [~, NEvents]=size(DataStruct.Events);
    isBound = 1;
    for BoundPos=1:NEvents
        if strcmp(user_event,DataStruct.Events(BoundPos).label)
            break
        elseif BoundPos==NEvents
            strMsg = ['No ' user_event ' event found in the data, processing all data points.'];
            bst_report('Warning', sProcess, [], strMsg);
            isBound = 0;
        end
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
   
 
   
  [itc, interval, Nwin]=Compute(sProcess,sInput,inputData, WindowLength, imagingKernel);
  
   %---Power Spectrum---
   itc = permute(itc, [1 3 2]); 
   %Output file creation
   f=interval*sRate/WindowLength;
   [ChannelNumber,~ ,~]=size(itc);
   Rows=1:1:ChannelNumber;
   FileMat = db_template('timefreqmat');
   OutputFiles=bst_process('GetNewFilename',bst_fileparts(Filepath),'timefreq_psd');
   FileMat.ChannelFlag=DataStruct.ChannelFlag;
   FileMat.TF=itc;
   FileMat.Comment= sprintf('ITC: %d/%d ms %s',sum(Nwin),WindowLength*(1/sRate)*1000,'Phase');
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
   
   FileMat.Measure='other';
   FileMat.Method='ITC';
   FileMat.DataFile=sInput.FileName;
   FileMat.nAvg=DataStruct.nAvg;
   FileMat.Options=sProcess.options;
   FileMat.History=DataStruct.History;
   FileMat.DisplayUnits = 'No Units';
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
   FileMat = bst_history('add', FileMat, 'compute', 'Frequency Tagging Analysis - Inter Trial Phase Coherence');
   % Save the new file
    save(OutputFiles,'-struct','FileMat');
    % Reference OutputFile in the database:
    db_add_data(sInput.iStudy, OutputFiles,FileMat);
end

%% ===== COMPUTE =====
function [itc, interval, Nwin]=Compute(sProcess,sInput,inputData, WindowLength, imagingKernel)
    
    [itc,interval,Nwin] = fta_itc(inputData,WindowLength,imagingKernel);
    bst_report('Info', sProcess, sInput, sprintf('Number of windows used: %d\n',Nwin));
    
end

%% ===== Computation logic for ITC goes here =====

function [itc,interval,N] = fta_itc(data,wl,ImagingKernel)
    % [itc,an,f] = fta_itc(data,wl,srate)
    % Computes Inter-Trial Coherence on all channels of EEG data over
    % consecutive non-overlapping windows (trials)
    %
    % Inputs:
    % data = cell structure with data{n}=EEG.data from EEGLAB dataset relative to epoch n
    % wl = length of the sliding window for power spectrum computation (in time points)
    % srate = data sampling rate
    %
    % Outputs:
    % itc = inter-trial coherence (channels x frequency x epochs)
    % interval = time interval points to compute the frequency vector
    %
    % Author: Marco Buiatti, CIMeC (University of Trento, Italy), 2016-2017.
    nep=length(data);
    N=0; % trial counter
    for ep=1:nep
        if size(data{ep},2) >= wl
            n_loc=floor(size(data{ep},2)/wl);
            for tr=1:n_loc
                data_loc=data{ep}(:,(1+(tr-1)*wl:tr*wl));
                % ImagingKernel here
                if ~isempty(ImagingKernel)
                    data_loc = ImagingKernel * data_loc;
                end
                an_loc=angle(fft(data_loc,[],2));
                if floor(wl/2)==wl/2
                    an(:,:,N+tr)=an_loc(:,1:wl/2+1);
                else
                    an(:,:,N+tr)=an_loc(:,1:(wl-1)/2+1);
                end
            end
            N=N+n_loc;
        end
        bst_progress('inc',  ceil(1/nep*100));
    end
    
    itc=abs(mean(exp(1i*an),3));


    if floor(wl/2)==wl/2
        interval=0:1:(wl/2);
    else
        interval=0:1:((wl-1)/2);
    end
    
end
