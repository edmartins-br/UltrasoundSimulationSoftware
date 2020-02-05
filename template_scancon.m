%	template_scancon.m
%	template script for converting from r-sin(theta) data to x-y image
%	must load r-sin(theta) data: rsdata

% --> QUESTION I <--
% Scan convert the r-sin(theta) buffer to produce a sector scan image.
% Use bilinear interpolation to compute the image values on the
% sector scan image grid.  Matlab's "interp2" function will help you
% do bilinear interpolation.

% compute values needed for interpolation
[R,S] = ndgrid(r,sin_theta);
xx = linspace(-30,30,512);
zz = linspace(0,60,512);
[Z,X] = ndgrid(zz,xx);
RI = (X.*X+Z.*Z).^0.5;
SI = X./RI;

% Create image w/ bilinear interpolation

image = interp2(S, R, abs(rsdata), SI, RI, 'bilinear');
t = find(isnan(image));
image(t) = zeros(size(t));
image2 = interp2(S, R, abs(rsdata2), SI, RI, 'bilinear');
t = find(isnan(image2));
image2(t) = zeros(size(t));
image3 = interp2(S, R, abs(rsdata3), SI, RI, 'bilinear');
t = find(isnan(image3));
image3(t) = zeros(size(t));

% --> QUESTION I <--
% Use two images on a logarithmic scale to answer this question:
% one on a 40dB scale, the other on a 20dB scale

showimage(image, 0, 20); % Display 20 dB scale image
showimage(image, 0, 40); % Display 40 dB scale image
showimage(image2, 0, 20); % Display 20 dB scale image
showimage(image2, 0, 40); % Display 40 dB scale image
showimage(image3, 0, 20); % Display 20 dB scale image
showimage(image3, 0, 40); % Display 40 dB scale image
%
% ALso, to answer this question, go back and look at the
% image PRIOR to scan conversion, especially in contrasting
% Artifacts for the Full Aperture and Decimated Aperture
%
