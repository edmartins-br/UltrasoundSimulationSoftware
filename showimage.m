 

function hfigout = showimage(X, nrow, dBrange, newwin, gamma, nmap)
%function hfigout = showimage(X, nrow, dBrange, newwin, gamma, nmap)
%
%	display image X in grayscale
%
%	nrow 	if 0, display image in natural coordinates
%		if nonzero, use as number of rows to 'panelize' skinny images
%		if negative, use as a number of columns
%
%	dBmin	if 0, use traditional linear scale
%		if nonzero, use as range in deciBels for display

%	10/5/96	Fessler, based on code by M. Lubinski

% Check inputs
if nargin < 1
	help showimage
	return
end
if nargin < 2
	nrow = 0;
end
if nargin < 3
	dBrange = 0;
end
if nargin < 4
	newwin = 1;		% put in new figure window ?
end
if nargin < 5
	gamma = 1.8;		% Gamma correction for colormap
end
if nargin < 6
	nmap = 256;		% size of colormap
end

	% Scale input to colormap
	if dBrange == 0
		X = intscale(X, 1, nmap);
	else
		X = intscale(X, 1, nmap, 0, dBrange, 'db');
	end

	% convert to panel format
	if nrow ~= 0
		if nrow < 0
			nrow = ceil(size(X,1) / -nrow);
		end
		X = panelize(X, nrow, 1);
	end

	% create figure
	[nr,nc] = size(X);
	if newwin
		hfig = figure;
	else
		clf
		hfig = gcf;
	end

	% colormap
	cmap = gray(nmap).^(1/gamma);
%	cmap = 1-gray(nmap);
%	cmap = 1-gray(nmap).^(1/gamma);
%	cmap = (1-gray(nmap)).^(1/gamma);
	set(hfig, 'colormap', cmap)

	% display
%	haxis = axes('units', 'pixels', 'position', [1 1 nc nr]);
	haxis = axes('units', 'pixels');

	image(X)

	% Display axis off
	set(haxis, 'units', 'normalized', 'visible', 'off');
    axis('image')

%	set(haxis, 'aspectratio', [nc/nr 1])

	if nargout
		hfigout = hfig;
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function panel = panelize(X, nrow, addstrip)
%function panel = panelize(X, nrow, addstrip)
%
%	Rearrange a skinny matrix X (that is otherwise awkward to display)
%	by snipping columns and laying-out next to each other.
%	This will wrap the columns of X so arrays with many rows
%	can be displayed on the screen.
%
%	Note: panel is such that element X(1,1) is in the upper left
%	and wrapped columns appear in order going down the image

%	10/6/96	Fessler, based on code by Lubinski

% Check inputs
if nargin < 1
	help panelize
	error 'arguments'
end
if nargin < 2
	nrow = ceil(sqrt(prod(size(X))));
end
if nargin < 3
	addstrip = 1;		% add bright strip between columns
end

nrow = min(nrow,size(X,1));

% Bright strip to separate columns
if addstrip & size(X,1) > nrow
	xmax = max(X(:));
	X = [X xmax(ones(size(X,1),1))];
end

% Zeros at end to round out
if rem(size(X,1),nrow)
	if addstrip
		X = [X; xmax(ones(1,size(X,2)))];	% one bright row strip 
	end
	if rem(size(X,1),nrow)
		X = [X; zeros(nrow-rem(size(X,1),nrow),size(X,2))];
	end
end

% Put in multiple columns
	[nr,nc] = size(X);
	panel = zeros(nrow, nr*nc/nrow);
	for ic = [1:(nr/nrow)]-1
		panel(:,ic*nc+[1:nc]) = X(ic*nrow+[1:nrow],:);
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function idx = intscale(a, imin, imax, amin, amax, method)
%SCALE	Scales data into integers.
%	intscale(A, imin, imax)
%		scales the data in A to the range [imin,imax]
%
%	intscale(A, imin, imax, amin, amax)
%		scales after truncating the data in A to range [amin,amax]
%
%	intscale(A, imin, imax, amin, amax, method)
%		'linear' - linear scale
%		'db'	 - log scale in deciBels

%       [No Guarantees. M. Lubinski]
%	10/6/96	Fessler, imin, imax

% Check function call
if (nargin < 3) | (nargin == 4)
	help inscale
	error('Not enough input arguments')
end
if nargin < 6
	method = 'linear';
end

if ~isreal(a)
	disp('Warning: using magnitude of complex array')
	a = abs(a);
end

if strcmp(method, 'db')
	% set nonpositives to min positive
	t = find(a <= 0);
	if t
		amin = max(amin, min(a(a > 0)));
		a(t) = amin * ones(size(t));
	end
	a = 20 * log(a) / log(10);
	a = a - max(a(:));

	if nargin < 5
		amin = -20;
	else
		amin = -amax;
	end
	amax = 0;

elseif nargin < 5
	amin = min(a(:));
	amax = max(a(:));
end

if amin == amax
	error('min = max')
end

% Scale
a = min(a, amax);
a = max(a, amin);
idx = imin + round( (a-amin) / (amax-amin) * (imax-imin) );
