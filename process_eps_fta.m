function varargout = process_eps_fta(varargin )
%%function [TF, FreqVector, Nwin, Messages] = process_psd_fta_test( F, sfreq, WinLength, WinOverlap, BadSegments, ImagingKernel, isVariance, FileName, EventType )
eval(macro_method);
end
%% ===== GET DESCRIPTION =====
function sProcess = GetDescription()
    % Description the process
    sProcess.Comment     = 'Compute Power Spectrum';
    sProcess.FileTag     = 'fta_ps';
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
    % Options: condition
    sProcess.options.condition.Comment = 'Condition';
    sProcess.options.condition.Type    = 'text';
    sProcess.options.condition.Value   = 'Power';
    % Options: Boundary event Label
    sProcess.options.BL.Comment = 'Boundary Event Label';
    sProcess.options.BL.Type    = 'text';
    sProcess.options.BL.Value   = 'boundary';
    sProcess.options.BL.Hidden  = 1; % Not required
    % === Time window
    sProcess.options.windowlen.Comment = 'Window Length';
    sProcess.options.windowlen.Type    = 'value';
    sProcess.options.windowlen.Value   = {10, 'seconds', 3};
    % === Time points duration
    sProcess.options.tpduration.Comment = 'Time point duration';
    sProcess.options.tpduration.Type    = 'value';
    sProcess.options.tpduration.Value   = {0.004, 'seconds', 3};
    sProcess.options.tpduration.Hidden  = 1; % Not required
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
    sProcess.options.isZeroPad.Class       = 'cZeroPad';
    
    % === Zero-Padding
    sProcess.options.minlengthzeropad.Comment     = 'Minimum Epoch Length';
    sProcess.options.minlengthzeropad.Type        = 'value';
    sProcess.options.minlengthzeropad.Value       = {5, 'seconds', 1};
    sProcess.options.minlengthzeropad.Class       = 'cZeroPad';
    sProcess.options.minlengthzeropad.Hidden      = 1;
    
    
    
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

    % Compute Sampling Rate from the Time vector
    if (~isfield(DataStruct, 'Time') || isempty(DataStruct.Time))
       bst_report('Error', sProcess, sInput, 'No Time Points exist in the data');
       return;
    else
       dT  = mean(diff(DataStruct.Time));
       sRate  = round(1/dT);
       strMsg = ['The sampling rate is ' num2str(sRate) ''];
       bst_report('Info',    sProcess, sInput, strMsg);
    end
    
    %Window is positive and not =0
   WindowLength=sProcess.options.windowlen.Value{1}*sRate; % in samples
    
   % Zeropadding
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
   
   
   %Window < Datasize

   [ChannelNumber, Datasize]= size(DataStruct.F);
   if(isZeroPad && Datasize < WindowLength && Datasize >= min_length_zp{1}/Srate && Datasize <= WindowLength)
       %Perform zero-padding
       DataStruct.F = [DataStruct.F zeros(ChannelNumber, WindowLength - Datasize)];
   elseif (WindowLength>Datasize)
        bst_report('Error', sProcess, [], 'Too short data. Try Zero-Padding.');
        return;
   end
   
   psdKernel   = ones(1,WindowLength)/WindowLength; % Taper is Square
   
   
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
  
 
  [ps, interval, Nwin]=Compute(sProcess,DataStruct.F, WindowLength, UpperBound, LowerBound, psdKernel);
    

   
   
   %---Power Spectrum---
   ps = permute(ps, [1 3 2]); % why?
   %Output file creation
   f=interval*sRate/WindowLength;
   [ChannelNumber,~ ,~]=size(ps);
   Rows=1:1:ChannelNumber;
   FileMat = db_template('timefreqmat');
   OutputFiles=bst_process('GetNewFilename',bst_fileparts(Filepath),'timefreq_psd');
   FileMat.ChannelFlag=DataStruct.ChannelFlag;
   FileMat.TF=ps;
   FileMat.Comment= sprintf('PSD: %d/%d s %s',Nwin,WindowLength*(1/sRate),sProcess.options.condition.Value);
   FileMat.DataType='data';
   FileMat.Time=[0,Datasize*(1/sRate)];
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
function [ps, interval, Nwin]=Compute(sProcess,subDataMat, WindowLength, UpperBound, LowerBound, psdKernal)
    
    [ChNumber, DataLength]= size(subDataMat);
    ps = zeros(ChNumber, WindowLength/2 +1); 
    
    w = WindowLength*(sum(psdKernal.^2));	%window squared and summed
    
    for ch = 1:size(subDataMat,1)
        if DataLength-WindowLength<floor(WindowLength*(1-UpperBound))  % use one window only if residual segment is shorter than wl/4
            DataLength=WindowLength;
            Nmax=ceil(2*DataLength/(WindowLength*LowerBound)); % max number of consecutive HALF windows (non-overlapping if l is multiple of wl/2)
            winStart=0;
        else
            Nmax=ceil(DataLength/(WindowLength*LowerBound)); % max number of consecutive HALF windows (non-overlapping if l is multiple of wl/2)
            winStart=floor((DataLength-WindowLength)/(Nmax-2)); % step
        end
        
        Nwin=Nmax-1; % number of consecutive full windows
        %disp(Nwin);
        ap_loc=0; %average power
        for i=1:Nwin
            wd=psdKernal.*subDataMat(ch,(i-1)*winStart+1:(i-1)*winStart+WindowLength);
            ap_loc=ap_loc + abs(fft(wd)).^2;
        end
        
        ap = ap_loc; % output average power
        
        n=sum(Nwin);
        ps(ch,1)= ap(1)/(w*n);
        ps(ch,2:WindowLength/2)=(ap(2:WindowLength/2)+ap(WindowLength:-1:WindowLength/2 +2))/(w*n);
        ps(ch,WindowLength/2 +1)=ap(WindowLength/2 +1)/(w*n);
        interval=0:1:(WindowLength/2);   
    end
    
    
     bst_report('Info', sProcess, [], sprintf('Number of windows used: %d',Nwin));
end
