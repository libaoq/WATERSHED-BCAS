function [widths, angles] = octopusDoppler_step1(RRslice, ilst, indices)

%ilst = 40:140; % z indices for motion correction
%indices =  60:100; % these are the Z indices to calculate doppler on

% I am really not sure what these are
step = 2;      % I am assuming this is the number of A scans per lateral resolution 
sz = 4 * step; % this is close to right, copied from getNAscan.m I don't know what it is

%% motion correction from octopus code
imgc = RRslice;
for i=1:1:size(imgc,2)-1
    r = corrcoef(imgc(ilst,i),imgc(ilst,i+1)); % add specific indices here.
    ccmat(i) = r(1,2);                         % ccmat is unused, Baoqiang 
    bulk_motion_estimate(i) = angle(r(1,2));
end
phase_diff_vect1 = [0;medfilt2(bulk_motion_estimate,[1,1],'symmetric')']; % phase difference calculated directly. "Exact" motion correction
phase_diff_vect2 = [0;repmat(0.0,[size(imgc,2)-1,1])];                    % phase difference assuming linear bulk motion. but not used, Baoqiang

phase_shift_vect = cumsum(phase_diff_vect1); % phase shift relative to the first line.
bulk_motion_corrector = exp(1i*phase_shift_vect)';
imgc2  = imgc.*repmat(bulk_motion_corrector,[size(imgc,1),1]);

%% high-pass filtering? take the diff between adjacent A-lines, if not at the exact same location, static material will not be efficiently diminished 
%  has to evaluate its impact on flow estimation
% h = [1,-1];
% imgc3 = conv2(imgc2,h,'same'); % used to remove stationary material. should be possible to improve filtering.

% high-pass filtering
tacq = 1/46780; 
h = design_filter_kernel(4*tacq,tacq); % high-pass filter taken from Fred's code
imgc3 = conv2(imgc2,h,'same'); 

%% other steps from octopus code, haven't yet figured out what they are for
maxlag = 2*sz+1;
    
transverse_indices = 1:size(imgc3,2); % these are the transverse indices to calculate on

img2xcorr = imgc3;

% in case there are reflections take noise variance at two locations.
yy = size(imgc3,2); 
noisevar1 = mean(var(imgc3(indices,1:round(yy*0.25))))*(2*sz+1)/2; % don't understand why multiply by (2*sz+1)/2
noisevar2 = mean(var(imgc3(indices,round(yy*0.75):yy)))*(2*sz+1)/2;
noisevar = min(noisevar1,noisevar2);

widths = zeros(length(indices),length(transverse_indices));
angles = zeros(length(indices),length(transverse_indices));
xcorrimg = zeros(length(transverse_indices),2*maxlag+1);

% hwait = waitbar( 0, sprintf('Calculating Z index 0 of %d',length(indices)) );
for i=1:1:length(indices)
    % waitbar(i/length(indices),hwait,sprintf('Calculating Z index %d of %d',i,length(indices))); 
    
    img_nolag = img2xcorr(indices(i),:);
    for ii=1:1:(2*maxlag+1)  % loop to create xcorr matrix for given axial position 
        shift = ii-maxlag-1; % shift goes from -maxlag to maxlag
        
        if shift > 0
            img2xcorr_shift = [zeros(1,shift),img2xcorr(indices(i),1:(length(img2xcorr)-shift))];
        else
            img2xcorr_shift = [img2xcorr(indices(i),(-shift+1):length(img2xcorr)),zeros(1,-shift)];
        end
        
        img_lag = conj(img2xcorr_shift);
        tempimg = conv_w_trunc((img_nolag.*img_lag)',ones(maxlag,1))';           % summation
        xcorrimg(:,ii) = tempimg(transverse_indices)*(maxlag-abs(shift))/maxlag; % why multiply by (maxlag-abs(shift))/maxlag
    end
    
    for ii=1:1:length(transverse_indices)
        min_index = transverse_indices(ii)-sz;
        max_index = transverse_indices(ii)+sz;
        
        if (min_index >= 1 && max_index <= size(img2xcorr,2))
            xcorrimg(ii,maxlag+1) = xcorrimg(ii,maxlag+1) - 2*noisevar;
            xcorrimg(ii,maxlag+2) = xcorrimg(ii,maxlag+2) + 1*noisevar;
            xcorrimg(ii,maxlag)   = xcorrimg(ii,maxlag)   + 1*noisevar;
            % xcorrimg(ii,:) = xcorrimg(ii,:)/max(xcorrimg(ii,:)); %do not
            % uncomment this line!!!
        end
        
        if abs(xcorrimg(ii,maxlag+1)) > 0.75*noisevar
            widths(i,ii) = calc_width(xcorrimg(ii,:)',1)';
            z(i,ii) = abs(xcorrimg(ii,maxlag+1)) < abs(xcorrimg(ii,maxlag)); %test if center value is less than next value
        end
    end
    
    %indicesmat = repmat(1:1:size(xcorrimg,2),[size(xcorrimg,1),1]);
    
    %OLD METHOD
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ii=1:1:size(xcorrimg,1)
        goodind = maxlag+1-ceil(widths(i,ii)/2) : maxlag+1+ceil(widths(i,ii)/2);
        %if ~isempty(goodind) & length(goodind) ~= max(goodind)-min(goodind)+1 error('ewr'); end
        if abs(xcorrimg(ii,maxlag+1)) > 0.75*noisevar
            x = (goodind-(maxlag+1))';
            %beta = unwrap(angle(xcorrimg(ii,goodind)))';
            %[Y,I] = max(abs(fftshift(fft(xcorrimg(ii,goodind),256,2),2)));
            angles(i,ii) = angle(xcorrimg(ii,maxlag+2)); %(1/(x'*x)*x'*beta);%(I-2^7-1) * pi/(2^7);%%
        else
            % angles(i,ii) = NaN;
            x = (goodind-(maxlag+1))';
            angles(i,ii) = angle(xcorrimg(ii,maxlag+2));
        end
        %widths2(i,:) = sqrt(sum(abs(xcorrimg).*indicesmat.^2,2)./sum(abs(xcorrimg),2) - (maxlag + 1)^2)';
    end
    
    %%%%%%%%%%%%%%%%
    %    xcorrimgfft = fftshift(fft(xcorrimg,256,2),2);%%%%%%%
    
    %      for ii=1:1:size(xcorrimg,1)
    %             if (xcorrimg(ii,maxlag+1)) > 0.75*noisevar
    %                 [Y,I]= max(abs(xcorrimgfft(ii,:)));%%%%CHANGE
    %                 angles(i,ii) = I - 2^7-1;
    %             else
    %                 angles(i,ii) = 0;
    %             end
    %      end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%x
end

%% phase unwrapping
ifunwrap = 0;
if ifunwrap == 1
    doppler_angle_unwrapped = unwrap(angles);
    difference = sum(doppler_angle_unwrapped~=angles);
    positions_to_unwrap = find(difference~=0);

    % This is a more careful way of unwrapping that makes sure the operation
    % is done in both directions, from Fred's code
    for i = positions_to_unwrap
        angles(:,i) = unwrap_single_line(angles(:,i));
    end  
end

% close(hwait)
