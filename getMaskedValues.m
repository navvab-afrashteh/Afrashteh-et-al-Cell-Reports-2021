function values = getMaskedValues (data,mask)
[r, c] = find(mask==1);
maskI = sub2ind(size(mask),r,c);
values = zeros(size(data,3),length(maskI));
for ii = 1:size(data,3)
    temp = data(:,:,ii);
    values(ii,:) = transpose(temp(maskI));
end