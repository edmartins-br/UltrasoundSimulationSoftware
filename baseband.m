 function bb = baseband(rf, df0)
%function bb = baseband(rf, df0)
%	Compute the baseband of RF data.
%	bb = BASEBAND(rf,df0) where df0 is digital center freq (0-0.5)
%	i.e. df0 = f0/fs = center_frequency / sampling_frequency
%
%	If RF is a matrix, baseband(RF,df0) does each column independently.
%	10/5/96	Fessler	simplified from Lubinski, added zero-padding

if isstr(rf) & strcmp(rf, 'test')
	test_baseband
	return
end

chat = 1;

[norig,M] = size(rf);
npad = norig;
if rem(log(norig)/log(2),1)
	npad = 2 ^ ceil(log(norig)/log(2));
	disp(sprintf('padding %g to length %g', norig, npad))
end

% Determine baseband shift index k0
k0 = round(npad * df0) + 1;
if chat
	disp(sprintf('k0 = %d, ideal=%g', k0, round(npad*df0)+1))
end

if (k0 > npad/2)
	disp(sprintf('k0 was %d, converting neg. freq. to pos. freq.',k0))
	k0 = npad-k0+1;
end

%
%	DFT
%
f = fft(rf, npad);

%
%	shift fft by k0 (& filter out neg. freq.)
%
b = zeros(npad,M);
%b((npad/2+2-k0):(npad+2-k0),:) = f(1:(npad/2+1),:);	% wrong sign!
b([npad/2+1:end]+1-k0,:) = f(npad/2+1:end,:);

% Shift fft back to flip halves of fft
for i = 1:M
	b(:,i) = fftshift(b(:,i));
end

% Take IDFT, multiplying by 2 to conserve power of signal
bb = ifft(b) * 2;
bb = bb(1:norig,:);

