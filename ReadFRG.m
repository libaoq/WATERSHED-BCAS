% Read FRG files
%{
filePath = '~/discostu4/R003/_raw/0531/053101.FRG';
frameNoAry = 1:2;
%}

function [FRG, ET] = ReadFRG(filePath, overscan, frameNoAry)		% [nk nx nframe]

% overscan always 1, suggested by thorlabs engineer, don't know what does it mean.
if nargin < 2
	overscan = 1; 
end

% Read header of the .FRG data using file handle fid
fid = fopen(filePath, 'rb' );
file_id = fread(fid, 16/1, 'uchar');        % File identification string
file_head = fread(fid, 24/4, 'int32');      % File_header
junk = fread(fid, 472/4, 'int32');             % reserved bytes

nfr = file_head(1);
nx = file_head(2);
nz = file_head(3);
% ny = file_head(4);		% #y frame
% nt = file_head(5);		% #vol
nk = file_head(6);
if nk == 0
    nk = nz*2;
end

% frameNoAry is used to choose the frames to reconstruct, e.g. [1 3 5 ...]
% to reconstruct only the odd frame
if nargin < 3
    frameNoAry = 1:nfr;
elseif frameNoAry(end) > nfr
    frameNoAry = frameNoAry(1):nfr;
end
		
% Read FRG
FRG = zeros(nk,nx,length(frameNoAry),'single');
ET = zeros(1,length(frameNoAry));

jfr = 0;
for ifr=1:nfr
	if length(find(frameNoAry == ifr)) > 0
        elapsed_time = fread(fid, 4/4, 'int32') /1000;	% msec
        frame_info = fread(fid, 36, 'int8');			% system time & reserved
        frame_data = fread(fid, nk*nx, 'int16');
        jfr = jfr+1;
        FRG(:,:,jfr) = reshape(frame_data, [nk nx]);
        ET(jfr) = elapsed_time;
        if (mod(jfr,ceil(length(frameNoAry)/2)) == 0)  
            disp(['... ReadFRG ' num2str(jfr) '/' num2str(length(frameNoAry)) '	' datestr(now,'HH:MM')]);  
        end
    else
        fseek(fid, 1*4+36+nk*nx*2, 'cof');
    end

    if ifr > max(frameNoAry)
        break;
    end
end

fclose(fid);

if overscan == 1
    % FRG = FRG(29:end,:,:); uncommented by Baoqiang
end

