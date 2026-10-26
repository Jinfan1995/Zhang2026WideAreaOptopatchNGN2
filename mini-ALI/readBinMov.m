
function [mov, nframe] = readBinMov(fileName, nrow, ncol, binDepth)
% [mov nframe] = readBinMov(fileName, nrow, ncol)
% If only one input: mov = readBinMov(fileName, nrow, ncol)
% The output, mov, is a 3-D array or unsigned 16-bit integers
% The binary input file is assumed to be written by Labview from a 
% Hamamatsu .dcimg file.  The binary file has 16-bit integers with a little
% endian format.
% The frame size information must be independently documented in a .txt 
% file and read into this function as the number of rows (nrow) and the 
% number of columns (ncol).
% Optional input to toggle 16bit (default) or 8bit inputs. -D.G. Itkis.

%set binDepth
if nargin < 4 || isempty(binDepth)
   binDepth = 16; 
end

%to avoid overflow when calculating total size of movie array.
nrow = uint64(nrow); 
ncol = uint64(ncol);

% read file into tmp vector
fProp = dir(fileName);
fSize = fProp.bytes;
fid = fopen(fileName);                  % open file

if binDepth == 16
    L = floor(fSize*(nrow*ncol));
    format = '*uint16';
elseif binDepth == 8
    L = floor(fSize*2*(nrow*ncol));
    format = '*uint8';
else 
    warning("Invalid input for bitDepth. Allowed values are 8 or 16. Defaulting to 16bit.");
    L = floor(fSize*(nrow*ncol));
    format = '*uint16';

end

mov = fread(fid,L*nrow*ncol, format, 'l');       % uint16, little endian
%mov = fread(fid,4000*nrow*ncol, format, 'l');  
fclose(fid);                            % close file

% reshape vector into appropriately oriented, 3D array
% L = int32(floor(length(mov)/(nrow*ncol)));
% tmp = tmp(1:L*nrow*ncol,1);
mov = reshape(mov,ncol, nrow, []);
mov = permute(mov, [2 1 3]);

if nargout > 1
    nframe = L;
end