function [ps,f,n]=fta_power_spectrum(data,wl,srate,overlap_factor,isPad,rejLength)
% fta_power_spectrum Computes power spectrum on all channels of EEG data.
%
% [ps, f, n] = fta_power_spectrum(data, winLength, sRate,isPad,padLength)
%
% Computes the power spectrum on all channels of EEG data by averaging over
% consecutive windows with overlapping, covering all data points.
%
% Inputs:
%   data  - Cell structure with data{n} = EEG data matrix (channels x time)
%           relative to epoch n (data{n} = EEG.data for EEGLAB datasets).
%   wl    - Length of the sliding window for power spectrum computation
%           (in time points).
%   sRate - Data sampling rate.
%   isPad - Turns on or off zero-padding option (Boolean = 0 or 1).
%   rejLength - If isPad is set to 1, an epoch with less than this amount
%               of data (in samples) will not be considered for zero-padding.
%               This denotes the minimum data length required to keep the
%               epoch.
%
% Outputs:
%   ps - Power spectrum (channels x frequency).
%   f  - Frequency bins.
%   n  - Number of epochs used to compute ps.
%
% Author: Marco Buiatti, CIMeC (University of Trento, Italy), 2016-.

% Adding zero-padding logic here

if(isPad)
    if(isempty(rejLength))
        rejLength = 1;
    end
    inputData_copy = data; % use the copy to perform zero-padding
    inputData_pad = [];
    count = 1; % Count of zero-padded (new) segments
    for iSeg=1:length(inputData_copy)
        [ChannelNumber, Datasize]= size(inputData_copy{iSeg});
        if(Datasize <= wl && Datasize >= rejLength)
            inputData_pad{count} = [inputData_copy{iSeg} zeros(ChannelNumber, wl - size(inputData_copy{iSeg},2))]; %Zero-padding
            strMsg = ['Data Segment : ' num2str(iSeg) ' is zero-padded'];
            disp(strMsg);
        else
            inputData_pad{count} = inputData_copy{iSeg};
        end
        count = count + 1;
    end
    data = inputData_pad;
end

% remove data epochs shorter than wl
for ep=1:length(data)
    loc(ep)=length(data{ep}(1,:))/wl;
end
data(loc<1)=[];

if isempty(data)
    error('Window longer than data!');
end

% Taper option (default taper is square)
option='s';
if  strcmp(option, 's')
    w=ones(wl,1)/wl;
elseif  strcmp(option, 'w')
    w=welch(wl);
elseif  strcmp(option, 'p')
    w=parzen(wl);
elseif  strcmp(option, 'h')
    w=hamming(wl);
end
w=w';

wss=wl*(sum(w.^2));	%window squared and summed

% compute ps for each electrode
overlap_length = round(wl * overlap_factor);

for el=1:size(data{1},1)
    ap=0;
    for ep=1:length(data)
        l=size(data{ep}(el,:),2);
        n_loc(ep) = floor((l - overlap_length) / (wl - overlap_length));
        
        ap_loc=0;
        % for each window, compute ps and sum to one-epoch ps
        for i=1:n_loc(ep)
            start_idx = (i - 1) * (wl - overlap_length) + 1;
            end_idx = start_idx + wl - 1;
            wd=w.*data{ep}(el,start_idx:end_idx);
            fwd = fft(wd); % FFT of wd 
            % fwd = fft(wd, [], 2); if the above line throws error, please comment that out and uncomment this line
            if ~isempty(imagingKernel)
                fwd = imagingKernel * fwd; % Computing the source PSD
            end
            ap_loc=ap_loc + abs(fwd).^2;
        end
        % sum to total ps
        ap=ap+ap_loc;
    end

    % normalize ps
    n=sum(n_loc);
    ps(el,1)=ap(1)/(wss*n);
    if floor(wl/2)==wl/2
        ps(el,2:wl/2)=2*ap(2:wl/2)/(wss*n);
        ps(el,wl/2 +1)=ap(wl/2 +1)/(wss*n);
    else
        ps(el,2:(wl+1)/2)=2*ap(2:(wl+1)/2)/(wss*n);
    end
    interval=0:1:(wl/2);
end
f=interval*srate/wl;

for ep=1:length(n_loc)
    disp(['Number of windows used for epoch ' num2str(ep) ' = ' num2str(n_loc(ep))]);
end