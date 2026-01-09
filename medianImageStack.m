function result = medianImageStack(source)

result=[];
X=size(source, 1);
Y=size(source, 2);
Z=size(source, 3);
temp=zeros(X+2, Y+2, Z+2);
temp(2:X+1, 2:Y+1, 2:Z+1)=source;

for i=2:Z+1
    temp(1, :, i)=temp(2, :, i);
    temp(X+2, :, i)=temp(X+1, :, i);
    temp(:, 1, i)=temp(:, 2, i);
    temp(:, Y+2, i)=temp(:, Y+1, i);    
end
temp(:,:,1)=temp(:,:,2);
temp(:,:,Z+2)=temp(:,:,Z+1);

% if ispc, % if Win32 computer, use dll file to process data
%     disp(['using dll file...']);
%     tt=median3d(uint16(temp));
% else     % if Unix or Win64, use slow procedure...
    disp(['using slow 3d filter file...']);
    tt = uint16(SavaMedian3D(temp));
% end;


result=tt(2:X+1, 2:Y+1, 2:Z+1);
