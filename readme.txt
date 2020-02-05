Files for ultrasound project

data06.mat	RF ultrasound data of Mystery Object (1500 x 65 matlab matrix)
template_beamform.m	A Matlab script template for doing the receive-only
			beamforming reconstruction (for questions A-H)
template_scancon.m	A Matlab script template for doing the scan conversion
              		(for questions I-K)
baseband.m	convert RF to baseband (see question E)
readme.txt	this file
showimage.m	image display utility, also converts skinny matrices to 'panels'

If you cannot copy the files, let me know ASAP and we'll work it out.

HELPFUL HINTS:
Useful Matlab functions. These functions may help you, but you don't 
necessarily HAVE to use them for this computer assignment.
     help     = get help on Matlab commands
     load     = load data into memory
     size     = get size of a matrix
     length   = get length of a vector
     ones     = create a vector or matrix filled with 1's
     zeros    = create a vector or matrix filled with 0's
     sqrt     = compute square roots
     exp      = compute exponentials
     abs      = compute absolute value (or magnitude of complex value)
     sum      = sum the elements of a vector/matrix
     angle    = compute phase of complex value
     find     = find indices of non-zero elements of a vector/matrix
     floor,ceil,round,fix = round or truncate values to integers
     linspace = created a linearly spaced vector
     meshgrid,ndgrid = transform 2 vectors into 2 array grids

If you have a column vector v of length n and you want to make a matrix M
of m columns all equal to v, i.e. M = [v v ... v], then the fastest way is:
	M = v(:,ones(1,m));

Warning! In Matlab :
       '  (single quote) is Hermitian transpose
       .' (period single quote) is transpose
      This is an important distinction when using complex values.

In Matlab, all math functions are geared toward matrix math, so
       *  (asterisk) is matrix multiplication
       .* (period asterisk) is element-by-element multiplication
      The same is true for / vs ./


Matlab runs fastest when the number of loops is minimized.
My solution has only one more loop (over transducer elements) beyond
the loop shown in the template.  While debugging you may want to reduce 
the number of beams and/or range
samples.

Good Luck everyone!
