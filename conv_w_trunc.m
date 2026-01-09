%This function convolves a matrix with a kerenel (2d), and truncates in the
%first dimension, based on the kernel length.  The returned matrix is the
%same size as the original matrix.
%USAGE: matret = conv_w_trunc(mat,kernel)
function matret = conv_w_trunc(mat,kernel)

matret = conv2(mat,kernel);
if size(kernel,1) >= 2 
    kernel_half_length = round(size(kernel,1)/2);
    matret = matret((kernel_half_length):(kernel_half_length -1+ size(mat,1)),:); %truncate to right number of pixels.
end
