function save_float_tiff(data, outdir, filename)
% SAVE_FLOAT_TIFF Save a signed floating-point array as a TIFF file.
%
%   save_float_tiff(data, outdir, filename)
%
%   data     : numeric array (single or double, can be signed)
%   outdir   : output directory (string or char)
%   filename : name of the file (without directory)
%
%   The TIFF will be saved as either 32-bit or 64-bit IEEE floating point.

    % Ensure directory exists
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end
    
    % Determine bit depth
    if isa(data, 'single')
        bitsPerSample = 32;
    elseif isa(data, 'double')
        bitsPerSample = 64;
    else
        % Convert everything else to single
        warning('Converting data to single precision for saving.');
        data = single(data);
        bitsPerSample = 32;
    end

    % Create full file path
    filepath = fullfile(outdir, filename);

    % Create Tiff object
    t = Tiff(filepath, 'w');

    % Set required tags
    tagstruct.ImageLength = size(data, 1);
    tagstruct.ImageWidth = size(data, 2);
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = bitsPerSample;
    tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP; % floating point
    tagstruct.SamplesPerPixel = 1;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Compression = Tiff.Compression.None;

    % Apply tags and write data
    t.setTag(tagstruct);
    t.write(data);
    t.close();
end
