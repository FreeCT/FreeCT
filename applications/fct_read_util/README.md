# FreeCT Reader Utility

This program interfaces with the FreeCT reader library to extract data from the TCIA LDCT dataset into a 2 file YAML+Binary format.  The two file approach is much faster to read for many applications.  We unfortunately do not have a fully working reader for the 2 file format, however that is in development.

**A enormous amount is in flux right now for FreeCT and I can't promise that everything in this document will be perfectly up to date (although I will do my best).  If you have questions or concerns, or notice something isn't working as it's supposed to, please reach out to johnmarianhoffman@gmail.com or create an issue in the issue tracker.  Thanks for your feedback and patience as I attempt to modernize FreeCT!**  - John 07/01/2020

# Building

Please follow the instructions for build FreeCT as a whole found in the README at https://github.com/FreeCT/FreeCT

# Running 

```bash
fct_read_util -i /path/to/dicom_directory/ -o /path/to/output_directory/
```

Note that while using absolute paths is always helpful, relative paths should work just fine!

### Input (-i)
The dicom directory should contain the TCIA LDCT projection data you wish to extract.  Other files or directories in the DICOM directory will likely cause the program to crash.  The good news is that as long as you have left the TCIA dataset as-is since it was downloaded, you should be able to extract the studies.

### Output (-o)

The program will create two files in the output directory:
 
```
/path/to/output_directory/meta.yaml
/path/to/output_directory/projections.dat
```
 
Metadata contains geometry and scan information in [YAML](yaml.org), a very convenient, human-readable, easily-edited serialization language.

The projection data is stored as a raw blob of single-precision floating point numbers.  Each projection's "header" data is 3 floats that store the source location, followed by n_detector_rows*n_detector_cols floats that store the actual projection data.

NOTE: In an update we will push soon, the header will be changed to store 3 floats for the source position, and 3 floats for any flying focal spot offset data.  We will update this README when that happens.  

A quick sample projection data reader for MATLAB might look something like:

```MATLAB
function [source_positions,projection_data] = readFreeCTProjectionData(filepath,
                                                                       start_projection,end_projection,
                                                                       n_detector_rows,n_detector_channels)
% filepath is /path/to/output_diretcory/projections.dat
% start_projection is zero-indexed projection
% end_projection is zero-indexed
% n_detector_rows for 3rd gen fanbeam CT is often 32, 64, or 128 (typically the smaller dimension)
% n_detector_channels is typically the larger detector dimension (>700 for Siemens, > 800 for GE)
% Note that MATLAB stores data in column-major arrays, while FreeCT assumes data is row-major.

sizeof_float     = 4; % bytes
frame_size_bytes = 3*sizeof_float + n_detector_rows*n_detector_channels*sizeof_float;
n_projections    = end_projection - start_projection + 1;

byte_offset      = start_projection * frame_size_bytes;
source_positions = zeros(n_projections,3);
projection_data  = zeros(n_detector_channels,n_detector_rows,n_projections);

fid = fopen(filepath,'r');
fseek(fid,byte_offset,'bof');
for (i = start_projection:end_projection)
  array_idx = i - start_projection + 1;
  source_positions(array_idx,:) = fread(fid,3,'float32');
  projection = fread(fid,n_detector_rows*n_detector_channels,'float32');
  projection_data(:,:,array_idx) = reshape(projection,n_detector_channels,n_detector_rows);
end

fclose(fid);

end
```
