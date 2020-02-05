%	template_beamform.m
%	Ultrasound beamforming project template
%
%	This script requires the variables
%	Data [ntime, nelem]	RF signals from each array element
%	f0			Transducer center frequency (MHz)
%	fs			sampling frequency (MHz)
%	c			speed of sound (mm/usec)
%	dx			transducer element spacing (mm)
%
% Get needed variables
%
load data06;
f0 = 4;		% MHz
fs = 16;	% MHz
c = 1.54;	% mm/us
% 
% Output Array Parameters
%
% Determine the array spacing dx in mm
lambda =c/f0;
dx = lambda/2;

deltat=1/fs;
[ntime, nelem] = size(Data);		% # time samples, # array elements
disp(sprintf('f0=%g MHz, deltat=%g usec, dx=%g mm', f0, deltat, dx))
disp(sprintf('# of Time Samples=%g,  # of Array Elements=%g',ntime,nelem))

%
% --> QUESTION A <--
% Make a wavefield plot of the raw Data
% Comment this out while you debug other parts of your program
%
showimage(Data, -4)
disp 'hit key', pause;
showimage(Data,-4,40)
disp 'hit key', pause;

%
% --> QUESTION B <--
%   Compute the number of total beams and beam spacing
%   I used variables called
%	nbeam = number of beams
%	sin_theta = vector of beam positions
%
nbeam = (((sqrt(2)/2)-(-(sqrt(2)/2))))/(1/65);
disp(sprintf('nbeam = %g', nbeam))
sin_theta = sin(1/65);
pause;

%
% It might be useful here to compute some values needed in the loop
%	below which will be constant inside the loop.
%	This will save cpu time in running your script.
%
?

%
% Create room for answers
%	This section merely makes matlab allocate memory before doing
%	any other processing. It is not a necessary section, but will
%	make Matlab run a little faster since it doesn't have to reallocate
%	space every time you add another beam
%
rsdata = zeros(ntime,nbeam);	% r-sin(theta) data buffer (part D)
rsdata2 = zeros(ntime,nbeam);	% r-sin(theta) data buffer (part G)
rsdata3 = zeros(ntime,nbeam);	% r-sin(theta) data buffer (part H)

% --> QUESTION E <--
% Convert Data to baseband
%
databb = baseband(Data, f0/fs);
%showimage(abs(databb), -4)
%disp 'hit key', pause

% --> QUESTION H <--
% Adjust the beam data according the windowing function
% databbwin = ??

%
% Do beamforming 
%
for ib=1:nbeam
	fprintf('Beam %d of %d\n', ib, nbeam)

	% --> QUESTION C,F <--
	% For the current beam, compute time delays (or delayed samples)
	% and phase rotations (QUESTION F) for each channel (i.e. transducer element)
	% as a function of range.
    
        for ie=1:nelem
                ??
        end

	% --> QUESTION D,G,H <--
	% Create data for r-sin(theta) buffer by computing the coherent
	% sum across the array using the delays (& phase rotations) from above
	?
	?
end

% Since rsdata is still at RF frequencies, we now convert it to baseband
rsdatabb = baseband(rsdata, f0/fs);
% Display results from part D 
% the contents of the r-sin(theta) buffer as 
% a gray scale image over a logarithmic scale of 40dB.
showimage(abs(rsdatabb), -4, 40)


% Display results from part G & H
showimage(abs(rsdata2), -4, 40)
showimage(abs(rsdata3), -4, 40)
