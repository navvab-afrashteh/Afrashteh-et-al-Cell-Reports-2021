function [h,str] = isSignificant(pvalue)
h = 1;
if pvalue<0.001
    str = '***';
elseif pvalue<0.01
    str = '**';
elseif pvalue<0.05
    str = '*';
else
    str = 'ns';
    h = 0;
end