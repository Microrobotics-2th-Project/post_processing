function [row,column,corr] = xcorr2_fftMine(Im, Tm)
%This function compute the cross-correlation between two matrices, searching
%for a template Tm into the given image Im, and return the (x,y)
%coordinates of the maximum (namely the center of the best match)
[Mt,Nt] = size(Tm);
[M,N] = size(Im);
If = fft2(Im);
Tf = fft2(Tm,M,N);
corr = abs(ifft2(If.*conj(Tf)));
[row,column] = find(corr==max(corr,[],'all')); 
row = mod(row+round(Mt/2),M);
column = mod(column+round(Nt/2),N);
if(row==0)
    row = M;
end
if(column==0)
    column = N;
end
end