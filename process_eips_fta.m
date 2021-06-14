function varargout = process_eips_fta(varargin )
%%function [TF, FreqVector, Nwin, Messages] = process_psd_fta_test( F, sfreq, WinLength, WinOverlap, BadSegments, ImagingKernel, isVariance, FileName, EventType )
eval(macro_method);
end
%% ===== GET DESCRIPTION =====
function sProcess = GetDescription()
    % Description the process
    sProcess.Comment     = 'Epoched Interlaced Power Spectrum';
    sProcess.FileTag     = 'ips';
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
    sProcess.options.sensortypes.Comment = 'Sensor types or names (empty=all): ';
    sProcess.options.sensortypes.Type    = 'text';
    sProcess.options.sensortypes.Value   = 'MEG, EEG';
    sProcess.options.sensortypes.InputTypes = {'data', 'results'};
    sProcess.options.sensortypes.Group   = 'input';
    % Options: condition
    sProcess.options.condition.Comment = 'condition';
    sProcess.options.condition.Type    = 'text';
    sProcess.options.condition.Value   = 'Power';
    % Options: Boundary event Label
    sProcess.options.BL.Comment = 'Boundary event label';
    sProcess.options.BL.Type    = 'text';
    sProcess.options.BL.Value   = 'boundary';
    % === Time window
    sProcess.options.epochduration.Comment = 'Epoch duration';
    sProcess.options.epochduration.Type    = 'value';
    sProcess.options.epochduration.Value   = {10, 'seconds', 3};
    % === Time points duration
    sProcess.options.tpduration.Comment = 'Time point duration';
    sProcess.options.tpduration.Type    = 'value';
    sProcess.options.tpduration.Value   = {0.004, 'seconds', 3};
    % === Upper bound overlap framing
    sProcess.options.Uepochoverlap.Comment = 'Upper bound epoch overlap';
    sProcess.options.Uepochoverlap.Type    = 'value';
    sProcess.options.Uepochoverlap.Value   = {75, '%', 0};
    % === Lower bound overlap framing
    sProcess.options.Lepochoverlap.Comment = 'Lower bound epoch overlap';
    sProcess.options.Lepochoverlap.Type    = 'value';
    sProcess.options.Lepochoverlap.Value   = {50, '%', 0};
    
    
end
%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) 
     
Comment = sProcess.Comment;

end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInput)
   % Initializing returned values
    OutputFiles = {};
    ImagingKernel={};
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
    ImagingKernel=RDataMat.ImagingKernel;
    Filepath=strtok(Filepath, '|');%path to the 'result' file, Rpath=path to the sensor file
    end
    %Singal rate (timepoint value) positive and not 0
    Srate=sProcess.options.tpduration.Value{1};
    if (Srate<=0)
       bst_report('Error', sProcess, [], 'Invalid timepoint value!');
        return;
   end
    %Window is positive and not =0
   WindowLength=2*floor((sProcess.options.epochduration.Value{1}/Srate)/2);%transformation from absolute time in datapoints of 0.004 seconds, needs to be an even integer
   if (WindowLength<1)
       bst_report('Error', sProcess, [], 'No time window selected!');
        return;
   end
   %Window < Datasize
   [ChannelNumber, Datasize]= size(DataStruct.F);
   if (WindowLength>Datasize)
       bst_report('Error', sProcess, [], 'Data shorter than Time window');
        return;
   end
   %lower limit<upper limit
   if (sProcess.options.Lepochoverlap.Value{1}>sProcess.options.Uepochoverlap.Value{1})
       bst_report('Warning', sProcess, sProcess.options.Lepochoverlap.Value{1}, 'Upper overlap bound lower than lower overlap bound, automatic switch');
        temp=sProcess.options.Lepochoverlap.Value{1};
        sProcess.options.Lepochoverlap.Value{1}=sProcess.options.Uepochoverlap.Value{1};
        sProcess.options.Uepochoverlap.Value{1}=temp;
   end
   %make the code easyer to read
  UpperBound=sProcess.options.Uepochoverlap.Value{1}/100;
  LowerBound=sProcess.options.Lepochoverlap.Value{1}/100;
  
   %---Extracting data Boundaries
   [~, NEvents]=size(DataStruct.Events);
   flag=0;
   for BoundPos=1:NEvents
       if strcmp(sProcess.options.BL.Value,DataStruct.Events(BoundPos).label)
       break
       elseif BoundPos==NEvents
           flag=1;
           bst_report('Warning', sProcess, [], 'No boundary event found, al data will be processed whitout epoch');
       end
   end
   
   %----executing if no Boundary found
    if flag
        [ap,TotWindowN]=Compute(sProcess,DataStruct.F, WindowLength, UpperBound, LowerBound, 0, ImagingKernel);
    else
        if isfield(DataStruct.Events, 'samples')
       Epochs=DataStruct.Events(BoundPos).samples;
        else
            Epochs=int32(DataStruct.Events(BoundPos).times/Srate);
        end
       [~,NEpochs]=size(Epochs);
       NEpochs=NEpochs-1;%the number of "epochs interval" is 1 less than the number of "epoch boundaries"
       %pre-allocating
       ap=zeros(ChannelNumber, WindowLength);	%average power
            %Apply ImagingKernel for right dimension of ap
            if ~isempty(ImagingKernel)
                ap = ImagingKernel * ap;
            end
       TotWindowN=0;
       
       for EpochCounter=1:NEpochs
           EpochedDataM=DataStruct.F(1:ChannelNumber, (Epochs(EpochCounter)+1):(Epochs(EpochCounter+1)-1));
           [ap_loc,WindowN_loc]=Compute(sProcess,EpochedDataM, WindowLength, UpperBound, LowerBound, EpochCounter, ImagingKernel);
           ap=ap+ap_loc;
           TotWindowN=TotWindowN+WindowN_loc;
       end
       
   end
   
   
   %---Power Spectrum---
   ps(:,1)=ap(:,1)/(TotWindowN);
   ps(:,2:WindowLength/2)=(ap(:,2:WindowLength/2)+ap(:,WindowLength:-1:WindowLength/2 +2))/(TotWindowN);
   ps(:,WindowLength/2 +1)=ap(:,WindowLength/2 +1)/(TotWindowN);
   ps = permute(ps, [1 3 2]);
   %Output file creation
   interval=0:1:(WindowLength/2);
   f=interval*(1/Srate)/WindowLength;
   [ChannelNumber,~ ,~]=size(ps);
   Rows=1:1:ChannelNumber;
   FileMat = db_template('timefreqmat');
   OutputFiles=bst_process('GetNewFilename',bst_fileparts(Filepath),'timefreq_psd');
   FileMat.ChannelFlag=DataStruct.ChannelFlag;
   FileMat.TF=ps;
   FileMat.Comment= sprintf('PSD: %d/%d ms %s',TotWindowN,WindowLength*Srate*1000,sProcess.options.condition.Value);
   FileMat.DataType='data';
   FileMat.Time=[0,Datasize*Srate];
   FileMat.Freqs=f;
   FileMat.RowNames=Rows;
   if(strcmp(sInput.FileType,'data')||strcmp(sInput.FileType,'raw'))
       FileMat.Device=DataStruct.Device;
       sRows={1,ChannelNumber};
       for ind=1:1:ChannelNumber
           sRows(1,ind)=cellstr(sprintf('E%d',ind));           
       end
      FileMat.RowNames=sRows;
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
   FileMat = bst_history('add', FileMat, 'compute', 'Interlaced Power Spectrum');
   % Save the new file
    save(OutputFiles,'-struct','FileMat');
    % Reference OutputFile in the database:
    db_add_data(sInput.iStudy, OutputFiles,FileMat);
end

%% ===== COMPUTE =====
function [ap,OWindowN]=Compute(sProcess,subDataMat, WindowLength, UpperBound, LowerBound, EpochCounter, ImagingKernel)
%inizializing values
ap=0;
OWindowN=0;
%----Determining Max and Min Step---

%Min step
   UMinStep=WindowLength*(1-UpperBound);%not rounded step size (brainstorm error in attempt to do it directly)
   MinStep=floor(UMinStep);
%Max step
   UMaxStep=WindowLength*(1-LowerBound);%not rounded step size (brainstorm error in attempt to do it directly)
   MaxStep=floor(UMaxStep);

%Dimension of current sub Data matrix
[ChannelNumber, Datasize]= size(subDataMat);
if (WindowLength>Datasize)
    error=sprintf('Epoch %d shorter than Time window',EpochCounter);
       bst_report('Warning', sProcess, [], error );
       return;
end
%Max number of windows
UMaxWindowN=((Datasize-WindowLength)/MinStep)+1; %not rounded number of windows (brainstorm error in attempt to do it directly)
MaxWindowN=floor(UMaxWindowN);

%Window number optimization
Gap=MaxStep-MinStep;
OStep=MinStep;
OWindowN=MaxWindowN;
Balancing=1;
if (UpperBound==LowerBound||Gap==0)
    error=sprintf('Epoch %d : No overlap gap, fixed overlap engaged, some data may not be considered.',EpochCounter);
    bst_report('Warning', sProcess, Gap, error);
    Balancing=0;
else
    while((Gap<(Datasize-(OStep*(OWindowN-1)+WindowLength)))&&OStep<MaxStep)
        OStep=OStep+1;
        UOWindowN=((Datasize-WindowLength)/OStep)+1; %not rounded number of windows (brainstorm error in attempt to do it directly)
        OWindowN=floor(UOWindowN);
    end
    if ((OStep==MaxStep)&&(Gap<(Datasize-(OStep*(OWindowN-1)+WindowLength))))
        error=sprintf('Epoch %d : Upper Bound too close from Lower Bound, some data is going to be left out as result.',EpochCounter);
        bst_report('Warning', sProcess, Gap, error);
        OStep=MinStep;
        OWindowN=MaxWindowN;
        Balancing=0;
    end
end

 %Pre-allocating
   WindowedDataM=zeros(ChannelNumber, WindowLength);
   ap=zeros(ChannelNumber, WindowLength);	%average power
   %Apply imaging kernel for right dimension of ap
    if ~isempty(ImagingKernel)
        ap = ImagingKernel * ap;
    end
   ChannelStart=1;
   ChannelEnd=WindowLength;
%%--START OF THE ACTUAL COMPUTING---
   for WindowNumber=1:(OWindowN-Balancing)       
       WindowedDataM=subDataMat(1:ChannelNumber, ChannelStart:ChannelEnd);
       WindowFFT=fft(WindowedDataM,[],2);
       %Apply imaging kernel
    if ~isempty(ImagingKernel)
        WindowFFT = ImagingKernel * WindowFFT;
    end
    ap=ap + abs(WindowFFT).^2;
       ChannelStart=ChannelStart+OStep;
       ChannelEnd=ChannelStart+WindowLength-1;
   end
   %Balancing
   if(Balancing==1)
       ChannelStart=(Datasize+1)-WindowLength;
       ChannelEnd=Datasize;
       WindowedDataM=subDataMat(1:ChannelNumber, ChannelStart:ChannelEnd);
       WindowFFT=fft(WindowedDataM,[],2);
       %Apply imaging kernel
    if ~isempty(ImagingKernel)
        WindowFFT = ImagingKernel * WindowFFT;
    end
    ap=ap + abs(WindowFFT).^2;
   end
      % 1 or more windows
   if (OWindowN==1)
       bst_report('Warning', sProcess, [], sprintf('Epoch Number %d: Data Only for one window',EpochCounter));
   end
end
