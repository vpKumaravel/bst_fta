function varargout = process_plfit_fta(varargin )
% process_plfit_fta: Given the power spectral density, the function computes the following:
%                       a) Normalized power spectrum @ Tagged-Frequency
%                       with respect to the user-defined frequency window
%                       b) Power spectrum @ Tagged-Frequency
%                       c) Power spectrum @ Background Frequency (or
%                       baseline)
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
% Authors: Marco Buiatti & Velu Prabhakar Kumaravel
%
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
    sProcess.options.NFR.Comment = ['Normalization Frequency Range ' char(177)];
    sProcess.options.NFR.Type    = 'value';
    sProcess.options.NFR.Value   = {0.3, 'Hz', 2};
    % === Tagged Frequency Alone
    sProcess.options.TFA.Comment = 'Power Spectrum @ Tagged Frequency';
    sProcess.options.TFA.Type    = 'checkbox';
    sProcess.options.TFA.Value   = 0;
    % === Baseline or Background Frequency Alone
    sProcess.options.BLA.Comment = 'Power Spectrum @ Baseline Frequencies';
    sProcess.options.BLA.Type    = 'checkbox';
    sProcess.options.BLA.Value   = 0;
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
    width=sProcess.options.NFR.Value{1};
    isTF = sProcess.options.TFA.Value;
    isBL = sProcess.options.BLA.Value;
    %check if TF is positive
    if(TaggedF<=0)
        bst_report('Error', sProcess, sInput, 'Selected Negative Tagged Frequency! Tagged Frequency must be positive');
        return;
    end
    %check if TF is less than the maximum frequency value
    if(TaggedF>DataStruct.Freqs(size(DataStruct.Freqs,2)))
        bst_report('Error', sProcess, sInput, 'Tagged Frequency greater than maximum data frequency');
        return;
    end
    %finding Tagged F position
    [~, PosTF]=min(abs(DataStruct.Freqs-TaggedF));   
    %finding upper and lower bound for frequency fit
    if (TaggedF-width)>0
    [~, PosLowBound]=min(abs(DataStruct.Freqs-(TaggedF-width)));
    else
        PosLowBound=2;%Because in position 1 we have freq 0 than in log=-inf and turns polyfit in NaN
    end
    
    if(TaggedF+width)<=max(DataStruct.Freqs)
        [~, PosUpBound]=min(abs(DataStruct.Freqs-(TaggedF+width)));
    else
        [~, PosUpBound]=size(DataStruct.Freqs, 2);
    end
    %check if frequency window is at least 3 frequency bin
    if(PosUpBound-PosLowBound<=3)
        bst_report('Warning', sProcess, sInput, 'Frequency window too small: Could lead to potentially unstable result');
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
            psbl(el)  =mean(DataMat(el,fitpoints));
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
        DataOut.HeadModelFile=DataStruct.HeadModelFile;
        DataOut.SurfaceFile=DataStruct.SurfaceFile;
        DataOut.nAvg=DataStruct.nAvg; 
        if(isfield(DataStruct, 'Whitener')); DataOut.Whitener=DataStruct.Whitener; end
        if(isfield(DataStruct, 'Function')); DataOut.Function=DataStruct.Function; end
        if(isfield(DataStruct, 'GoodChannel')); DataOut.GoodChannel=DataStruct.GoodChannel; end
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
    
    % Code for creating 2 other output files - one for tag frequency alone
    % and another for baseline or background frequencies alone
    if(isTF)
        %Creating output file for TF
        if strcmp(DataStruct.DataType,'data')
            DataOut= db_template('datamat');
            DataOut.F=zeros(size(DataMat,1),2);
            DataOut.F(:,1)=pstf;
            DataOut.F(:,2)=pstf;
            DataOut.Device=DataStruct.Device;
            OutputFiles=bst_process('GetNewFilename',bst_fileparts(sInput.FileName),'data_concat');
        elseif strcmp(DataStruct.DataType,'results')
            DataOut= db_template('resultsmat');
            DataOut.ImageGridAmp=zeros(size(DataMat,1),2);
            DataOut.ImageGridAmp(:,1)=pstf;
            DataOut.ImageGridAmp(:,2)=pstf;
            %extracting result path
            [ ~ ,Filepath]=strtok(sInput.DataFile, '|');%removing "link" text from string
            [Rpath,Filepath]=strtok(Filepath, '|');%storing the right filepath of the result
            RDataMat=in_bst_data(Rpath);
            DataOut.ImagingKernel=RDataMat.ImagingKernel;
            Filepath=strtok(Filepath, '|');%path to the 'result' file, Rpath=path to the sensor file
            OutputFiles=bst_process('GetNewFilename',bst_fileparts(Filepath),'results_concat');
            DataOut.DataFile=Filepath;
            DataOut.HeadModelFile=DataStruct.HeadModelFile;
            DataOut.SurfaceFile=DataStruct.SurfaceFile;
            DataOut.nAvg=DataStruct.nAvg;
            if(isfield(DataStruct, 'Whitener')); DataOut.Whitener=DataStruct.Whitener; end
            if(isfield(DataStruct, 'Function')); DataOut.Function=DataStruct.Function; end
            if(isfield(DataStruct, 'GoodChannel')); DataOut.GoodChannel=DataStruct.GoodChannel; end
        end
        
        % ===== SAVE FILE =====
        DataOut.ChannelFlag=DataStruct.ChannelFlag;
        DataOut.Comment=sprintf('TF %0.2f Hz',TaggedF);
        DataOut.Time=DataStruct.Time;
        DataOut.History=DataStruct.History;
        DataOut= bst_history('add', DataOut, 'compute', 'Frequency Tag Peak');
        DataOut.Options=sProcess.options;
        % Save the new file
        save(OutputFiles,'-struct','DataOut');
        % Reference OutputFile in the database:
        db_add_data(sInput.iStudy, OutputFiles,DataOut);
    end
    
    if(isBL)
        %Creating output file for BL
        if strcmp(DataStruct.DataType,'data')
            DataOut= db_template('datamat');
            DataOut.F=zeros(size(DataMat,1),2);
            DataOut.F(:,1) = psbl;
            DataOut.F(:,2) = psbl;
            DataOut.Device=DataStruct.Device;
            OutputFiles=bst_process('GetNewFilename',bst_fileparts(sInput.FileName),'data_concat');
        elseif strcmp(DataStruct.DataType,'results')
            DataOut= db_template('resultsmat');
            DataOut.ImageGridAmp=zeros(size(DataMat,1),2);
            DataOut.ImageGridAmp(:,1) = psbl;
            DataOut.ImageGridAmp(:,2) = psbl;
            %extracting result path
            [ ~ ,Filepath]=strtok(sInput.DataFile, '|');%removing "link" text from string
            [Rpath,Filepath]=strtok(Filepath, '|');%storing the right filepath of the result
            RDataMat=in_bst_data(Rpath);
            DataOut.ImagingKernel=RDataMat.ImagingKernel;
            Filepath=strtok(Filepath, '|');%path to the 'result' file, Rpath=path to the sensor file
            OutputFiles=bst_process('GetNewFilename',bst_fileparts(Filepath),'results_concat');
            DataOut.DataFile=Filepath;
            DataOut.HeadModelFile=DataStruct.HeadModelFile;
            DataOut.SurfaceFile=DataStruct.SurfaceFile;
            DataOut.nAvg=DataStruct.nAvg;
            if(isfield(DataStruct, 'Whitener')); DataOut.Whitener=DataStruct.Whitener; end
            if(isfield(DataStruct, 'Function')); DataOut.Function=DataStruct.Function; end
            if(isfield(DataStruct, 'GoodChannel')); DataOut.GoodChannel=DataStruct.GoodChannel; end
        end
        % ===== SAVE FILE =====
        DataOut.ChannelFlag=DataStruct.ChannelFlag;
        DataOut.Comment=sprintf('BL %0.2f Hz',TaggedF);
        DataOut.Time=DataStruct.Time;
        DataOut.History=DataStruct.History;
        DataOut= bst_history('add', DataOut, 'compute', 'Frequency Tag Peak');
        DataOut.Options=sProcess.options;
        % Save the new file
        save(OutputFiles,'-struct','DataOut');
        % Reference OutputFile in the database:
        db_add_data(sInput.iStudy, OutputFiles,DataOut);
    end
end