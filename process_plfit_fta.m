function varargout = process_plfit_fta(varargin )
%%function
eval(macro_method);
end
%% ===== GET DESCRIPTION =====
function sProcess = GetDescription()
    % Description the process
    sProcess.Comment     = 'Frequency Tag Peak';
    sProcess.FileTag     = 'plfit';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = 'Frequency Tagging Analysis';
    sProcess.Index       = 601;
    sProcess.Description = 'placeholder';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'timefreq'};
    sProcess.OutputTypes = {'data', 'results'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Options:
    % === Tagged Frequency
    sProcess.options.TaggedF.Comment = 'Tagged Frequency';
    sProcess.options.TaggedF.Type    = 'value';
    sProcess.options.TaggedF.Value   = {0.8, 'Hz', 2};
    % === Half Frequency Observation Window
    sProcess.options.HFOW.Comment = 'Half Frequency Observation Window';
    sProcess.options.HFOW.Type    = 'value';
    sProcess.options.HFOW.Value   = {1, 'Hz', 2};
    % === Power Law or Scalar Fit
    sProcess.options.CheckPowLaw.Comment = 'Power-law fit, if unchecked scalar fit will be applied';
    sProcess.options.CheckPowLaw.Type    = 'checkbox';
    sProcess.options.CheckPowLaw.Value   = 1;
    
    
end
%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)

Comment = sProcess.Comment;

end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInput)
    % Initialize returned values
    OutputFiles = {};
    %Loading Data Matrix
    DataStruct = in_bst_data(sInput.FileName);
    DataMat=squeeze(DataStruct.TF);
    %Inizialization
    np=ones(size(DataMat,1),1);
    TaggedF=sProcess.options.TaggedF.Value{1};    
    width=sProcess.options.HFOW.Value{1};
    %check if TF is positive
    if(TaggedF<=0)
    bst_report('Error', sProcess, [], 'Selected Negative Tagged Frequency! Tagged Frequency must be positive');
        return;
    end;
    %check if TF is less than the maximum frequency value
    if(TaggedF>DataStruct.Freqs(size(DataStruct.Freqs,2)))
    bst_report('Error', sProcess, [], 'Tagged Frequency greater than maximum data frequency');
        return;
    end;
    %finding Tagged F position
    [~, PosTF]=min(abs(DataStruct.Freqs-TaggedF))   
    %finding upper and lower bound for frequency fit
    if (TaggedF-width)>0
    [~, PosLowBound]=min(abs(DataStruct.Freqs-(TaggedF-width)));
    else
        PosLowBound=2;%Because in position 1 we have freq 0 than in log=-inf and turns polyfit in NaN
    end
    
    if(TaggedF+width)<=max(DataStruct.Freqs);
        [~, PosUpBound]=min(abs(DataStruct.Freqs-(TaggedF+width)));
    else
        [~, PosUpBound]=size(DataStruct.Freqs);
    end
    %check if frequency window is at least 3 frequency bin
    if(PosUpBound-PosLowBound<=3)
    bst_report('Error', sProcess, [], 'Frequency window too small!');
        return;
    end
    %determining frequency window
    fitpoints=[PosLowBound:PosTF-1 PosTF+1:PosUpBound];
    
    %computation
    fitorder=1;
    for el=1:size(DataMat,1)
        pstf(el)=DataMat(el,PosTF);
        if (sProcess.options.CheckPowLaw.Value)
            [Psp,Ssp] = polyfit(log(DataStruct.Freqs(fitpoints)),log(DataMat(el,fitpoints)),fitorder);
            psbl(el)  =exp(polyval(Psp,log(DataStruct.Freqs(PosTF))));
        else
            psbl(el)  =mean(DataMat(el,fitpoints))
        end
    end
    %Creating output file
    if strcmp(DataStruct.DataType,'data')
        DataOut= db_template('datamat');
        DataOut.F=zeros(size(DataMat,1),2);
        DataOut.F(:,1)=pstf./psbl;
        DataOut.F(:,2)=pstf./psbl;
        DataOut.Device=DataStruct.Device;
        OutputFiles=bst_process('GetNewFilename',bst_fileparts(sInput.FileName),'data_concat');
    elseif strcmp(DataStruct.DataType,'results')
        DataOut= db_template('resultsmat');
        DataOut.ImageGridAmp=zeros(size(DataMat,1),2);
        DataOut.ImageGridAmp(:,1)=pstf./psbl;
        DataOut.ImageGridAmp(:,2)=pstf./psbl;
        %extracting result path
        [ ~ ,Filepath]=strtok(sInput.DataFile, '|');%removing "link" text from string
        [Rpath,Filepath]=strtok(Filepath, '|');%storing the right filepath of the result
        RDataMat=in_bst_data(Rpath);
        DataOut.ImagingKernel=RDataMat.ImagingKernel;
        Filepath=strtok(Filepath, '|');%path to the 'result' file, Rpath=path to the sensor file
        OutputFiles=bst_process('GetNewFilename',bst_fileparts(Filepath),'results_concat');
        DataOut.DataFile=Filepath;
        DataOut.Function=DataStruct.Function;
        DataOut.HeadModelFile=DataStruct.HeadModelFile;
        DataOut.SurfaceFile=DataStruct.SurfaceFile;
        DataOut.nAvg=DataStruct.nAvg;
        DataOut.Whitener=DataStruct.Whitener;
        DataOut.GoodChannel=DataStruct.GoodChannel;
    end
    % ===== SAVE FILE =====  
    DataOut.ChannelFlag=DataStruct.ChannelFlag;
    DataOut.Comment=sprintf('NP %0.2f Hz',TaggedF);
    DataOut.Time=DataStruct.Time;
    DataOut.History=DataStruct.History;
    DataOut= bst_history('add', DataOut, 'compute', 'Frequency Tag Peak');
    DataOut.Options=sProcess.options;
    % Save the new file
    save(OutputFiles,'-struct','DataOut');
    % Reference OutputFile in the database:
    db_add_data(sInput.iStudy, OutputFiles,DataOut);
end
