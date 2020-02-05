function output=fwhm(in)
% finds fwhm of a smooth function
[y,i] = max(in);
j = find(in(1:i)>(y/2));
jj = j(1);
if jj > 1
  zc1 = jj - (in(jj)-y/2)/(in(jj)-in(jj-1));
else
  zc1 = 1;
end
j = find(in(i:length(in))<(y/2));
if length(j) < 1
  zc2 = length(in);
else
  jj = j(1)+i-1;
  zc2 = jj - 1 + (in(jj-1)-y/2)/(in(jj-1)-in(jj));
end
output = zc2-zc1;
