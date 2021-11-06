function B = tempFilt_imgStack(A,g)

[dim1, dim2, nFrames] = size(A);
if nFrames>1
    nPix = dim1*dim2;
    B = reshape(A,dim1*dim2,nFrames);
    parfor p = 1:nPix
        B(p,:) = filtfilt(g,1,double(B(p,:)));
    end
    B = reshape(B,dim1,dim2,nFrames);
else
    B = A;
end
