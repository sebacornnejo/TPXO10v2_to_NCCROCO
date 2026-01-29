
% Script to convert TPXO10v2 Atlas binary files to a NetCDF file 
% compatible with CROCO tools (following the TPXO7.nc structure).
%
% This script maintains the high resolution of TPXO10v2 (1/30 degree)
% but maps the variables to the names expected by legacy tools:
%   h (wct), ssh_r, ssh_i, u_r, u_i, v_r, v_i.
% It uses a two-pass approach to handle large memory requirements.
%
% ORIENTATION FIX:
% - XY returns S-N (increasing lat).
% - We do NOT flip lat or data. 
% - Result is strictly S-N (-90 to 90), compatible with CROCO tools.

clear all
close all
clc

%% Configuration
% Paths
tmd_path = './TMD2.5';
data_dir = '.tide_models/TPXO10_atlas_v2';
output_file = './DATASETS_CROCOTOOLS/TPXO7/TPXO10v2.nc';

con_string = 'M2 S2 N2 K2 K1 O1 P1 Q1 Mf Mm M4 2N2 Mn4 Ms4 S1';

% Add TMD path 
addpath(genpath(tmd_path));

% Temp directory for intermediate files
temp_dir = fileparts(output_file);
if isempty (temp_dir)
  temp_dir = pwd;
end

    % % 1. Read and Extend Grid disp('Reading Grid...');

% Read grid using TMD2.5 grd_in
[ll_lims, wct, mask, ~, ~] = grd_in(fullfile(data_dir, 'grid_tpxo10_atlas_30_v2'));

% Get coordinates
[lon, lat] = XY(ll_lims, size(wct, 1), size(wct, 2));
% lat = fliplr(lat);
% REMOVED : Keep S - N

                         % Get masks
[~, masku, maskv] = Muv(mask);

% Reorient and cast
% Transpose [Lon, Lat] -> [Lat, Lon]
% Do NOT flipud (keeps S-N)
mask = uint8(mask');
masku = uint8(masku');
maskv = uint8(maskv');
wct = double(wct'); 

% REMOVED EXTENSION: croco_pytools handles wrapping natively if grid is 0-360.
% Manually extending creates >360 span which breaks strict periodicity checks.
% lon = [lon(end)-360, lon, lon(1)+360];
% wct = cat(2, wct(:,end), wct, wct(:,1));
% mask = cat(2, mask(:,end), mask, mask(:,1));
% masku = cat(2, masku(:,end), masku, masku(:,1));
% maskv = cat(2, maskv(:,end), maskv, masku(:,1)); 

%% 2. Prepare Constituents
cons_cell = strsplit(con_string, ' ');
Ncons = length(cons_cell);

disp('Getting Constituent Parameters...');
try
    [ispec, amp, ph, omega, alpha] = tmd_constit(lower(cons_cell));
catch ME
    error('Error calling tmd_constit. Make sure the fixed tmd_constit.m is in your path. msg: %s', ME.message);
end

%% 3. Pass 1: Process and Save Temp
disp('Starting Pass 1: Processing Constituents...');

for k = 1:Ncons
    cname = cons_cell{k};
    cname_file = lower(cname);
    fprintf('  Processing %s (%d/%d)...\n', cname, k, Ncons);
    
    % --- H (Elevation) ---
    fname_h = fullfile(data_dir, ['h_', cname_file, '_tpxo10_atlas_30_v2']);
    tmp = h_in(fname_h, 1);
    
    % Transpose only, no flipud
    h_slice = tmp'; 
    
    % Extend - REMOVED
    % h_slice = cat(2, h_slice(:,end), h_slice, h_slice(:,1));
    
    % Regionfill NaNs/Land
    tmp_fill = complex(regionfill(real(h_slice), mask==0|isnan(h_slice)), ...
                       regionfill(imag(h_slice), mask==0|isnan(h_slice)));
    ssh_r_slice = real(tmp_fill);
    ssh_i_slice = imag(tmp_fill);
    
    % --- U/V (Transport) ---
    fname_u = fullfile(data_dir, ['u_', cname_file, '_tpxo10_atlas_30_v2']);
    [tmpu, tmpv] = u_in(fname_u, 1);
    
    % U
    u_slice = tmpu'; % No flipud
    % u_slice = cat(2, u_slice(:,end), u_slice, u_slice(:,1)); % REMOVED
    tmp_fill = complex(regionfill(real(u_slice), masku==0|isnan(u_slice)), ...
                       regionfill(imag(u_slice), masku==0|isnan(u_slice)));
    u_r_slice = real(tmp_fill);
    u_i_slice = imag(tmp_fill);
    
    % V
    v_slice = tmpv'; % No flipud
    % v_slice = cat(2, v_slice(:,end), v_slice, v_slice(:,1)); % REMOVED
    tmp_fill = complex(regionfill(real(v_slice), maskv==0|isnan(v_slice)), ...
                       regionfill(imag(v_slice), maskv==0|isnan(v_slice)));
    v_r_slice = real(tmp_fill);
    v_i_slice = imag(tmp_fill);
    
    % Save to temp file
    temp_file = fullfile(temp_dir, sprintf('temp_pass1_%d.mat', k));
    save(temp_file, 'ssh_r_slice', 'ssh_i_slice', ...
                    'u_r_slice', 'u_i_slice', ...
                    'v_r_slice', 'v_i_slice', '-v7.3');
    
    clear ssh_r_slice ssh_i_slice u_r_slice u_i_slice v_r_slice v_i_slice h_slice u_slice v_slice tmp tmpu tmpv tmp_fill;
end

%% 4. Initialize NetCDF
disp('Initializing NetCDF file...');

mode = netcdf.getConstant('NETCDF4'); 
ncid = netcdf.create(output_file, mode);

% Dimensions
dim_lon_r = netcdf.defDim(ncid, 'lon_r', length(lon));
dim_lat_r = netcdf.defDim(ncid, 'lat_r', length(lat));
dim_lon_u = netcdf.defDim(ncid, 'lon_u', length(lon));
dim_lat_u = netcdf.defDim(ncid, 'lat_u', length(lat));
dim_lon_v = netcdf.defDim(ncid, 'lon_v', length(lon));
dim_lat_v = netcdf.defDim(ncid, 'lat_v', length(lat));
dim_periods = netcdf.defDim(ncid, 'periods', Ncons);

% Coordinate Variables
% lon_r
var_lon_r = netcdf.defVar(ncid, 'lon_r', 'NC_FLOAT', dim_lon_r);
netcdf.putAtt(ncid, var_lon_r, 'long_name', 'Longitude at SSH points');
netcdf.putAtt(ncid, var_lon_r, 'units', 'degrees_east');

% lat_r
var_lat_r = netcdf.defVar(ncid, 'lat_r', 'NC_FLOAT', dim_lat_r);
netcdf.putAtt(ncid, var_lat_r, 'long_name', 'Latitude at SSH points');
netcdf.putAtt(ncid, var_lat_r, 'units', 'degrees_north');

% lon_u
var_lon_u = netcdf.defVar(ncid, 'lon_u', 'NC_FLOAT', dim_lon_u);
netcdf.putAtt(ncid, var_lon_u, 'long_name', 'Longitude at U points');
netcdf.putAtt(ncid, var_lon_u, 'units', 'degrees_east');

% lat_u
var_lat_u = netcdf.defVar(ncid, 'lat_u', 'NC_FLOAT', dim_lat_u);
netcdf.putAtt(ncid, var_lat_u, 'long_name', 'Latitude at U points');
netcdf.putAtt(ncid, var_lat_u, 'units', 'degrees_north');

% lon_v
var_lon_v = netcdf.defVar(ncid, 'lon_v', 'NC_FLOAT', dim_lon_v);
netcdf.putAtt(ncid, var_lon_v, 'long_name', 'Longitude at V points');
netcdf.putAtt(ncid, var_lon_v, 'units', 'degrees_east');

% lat_v
var_lat_v = netcdf.defVar(ncid, 'lat_v', 'NC_FLOAT', dim_lat_v);
netcdf.putAtt(ncid, var_lat_v, 'long_name', 'Latitude at V points');
netcdf.putAtt(ncid, var_lat_v, 'units', 'degrees_north');

% periods (optional, usually just index, but let's put indices)
var_periods = netcdf.defVar(ncid, 'periods', 'NC_FLOAT', dim_periods);
netcdf.putAtt(ncid, var_periods, 'long_name', 'Tide periods');

% Data Variables
% h (WCT)
var_h = netcdf.defVar(ncid, 'h', 'NC_FLOAT', [dim_lon_r dim_lat_r]);
netcdf.putAtt(ncid, var_h, 'long_name', 'Topography');
netcdf.putAtt(ncid, var_h, 'units', 'm');

% ssh_r
var_ssh_r = netcdf.defVar(ncid, 'ssh_r', 'NC_FLOAT', [dim_lon_r dim_lat_r dim_periods]);
netcdf.putAtt(ncid, var_ssh_r, 'long_name', 'Elevation real part');
netcdf.putAtt(ncid, var_ssh_r, 'units', 'm');

% ssh_i
var_ssh_i = netcdf.defVar(ncid, 'ssh_i', 'NC_FLOAT', [dim_lon_r dim_lat_r dim_periods]);
netcdf.putAtt(ncid, var_ssh_i, 'long_name', 'Elevation imaginary part');
netcdf.putAtt(ncid, var_ssh_i, 'units', 'm');

% u_r
var_u_r = netcdf.defVar(ncid, 'u_r', 'NC_FLOAT', [dim_lon_u dim_lat_u dim_periods]);
netcdf.putAtt(ncid, var_u_r, 'long_name', 'U-transport component real part');
netcdf.putAtt(ncid, var_u_r, 'units', 'm2.s-1');

% u_i
var_u_i = netcdf.defVar(ncid, 'u_i', 'NC_FLOAT', [dim_lon_u dim_lat_u dim_periods]);
netcdf.putAtt(ncid, var_u_i, 'long_name', 'U-transport component imaginary part');
netcdf.putAtt(ncid, var_u_i, 'units', 'm2.s-1');

% v_r
var_v_r = netcdf.defVar(ncid, 'v_r', 'NC_FLOAT', [dim_lon_v dim_lat_v dim_periods]);
netcdf.putAtt(ncid, var_v_r, 'long_name', 'V-transport component real part');
netcdf.putAtt(ncid, var_v_r, 'units', 'm2.s-1');

% v_i
var_v_i = netcdf.defVar(ncid, 'v_i', 'NC_FLOAT', [dim_lon_v dim_lat_v dim_periods]);
netcdf.putAtt(ncid, var_v_i, 'long_name', 'V-transport component imaginary part');
netcdf.putAtt(ncid, var_v_i, 'units', 'm2.s-1');

% Global Attributes
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'title', 'TPXO10v2 Atlas (Compat TPXO7)');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'date', datestr(now));
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'components', con_string);
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'description', 'Generated from TPXO10v2 Atlas binary files to match TPXO7 structure for CROCO tools. Orientation: S-N.');

netcdf.endDef(ncid);

%% 5. Write Data
disp('Writing Coordinates and Topography...');

netcdf.putVar(ncid, var_lon_r, lon);
netcdf.putVar(ncid, var_lat_r, lat);
netcdf.putVar(ncid, var_lon_u, lon);
netcdf.putVar(ncid, var_lat_u, lat);
netcdf.putVar(ncid, var_lon_v, lon);
netcdf.putVar(ncid, var_lat_v, lat);

% Calculate periods in hours
% omega is in radians/second
% Period (hours) = (2*pi / omega) / 3600
periods_hours = (2*pi ./ omega) / 3600;
netcdf.putVar(ncid, var_periods, periods_hours);

% Write wct (h)
netcdf.putVar(ncid, var_h, permute(wct, [2 1]));

disp('Starting Pass 2: Writing Constituents to NetCDF...');
for k = 1:Ncons
    fprintf('  Writing Constituent %d/%d...\n', k, Ncons);
    
    temp_file = fullfile(temp_dir, sprintf('temp_pass1_%d.mat', k));
    load(temp_file);
    
    start = [0 0 k-1];
    
    % Permute [Lat, Lon] -> [Lon, Lat] for NetCDF. 
    % Note: lat is dim 2, lon is dim 1 in matlab storage.
    count = [size(ssh_r_slice, 2), size(ssh_r_slice, 1), 1];
    
    netcdf.putVar(ncid, var_ssh_r, start, count, permute(ssh_r_slice, [2 1]));
    netcdf.putVar(ncid, var_ssh_i, start, count, permute(ssh_i_slice, [2 1]));
    
    netcdf.putVar(ncid, var_u_r, start, count, permute(u_r_slice, [2 1]));
    netcdf.putVar(ncid, var_u_i, start, count, permute(u_i_slice, [2 1]));
    
    netcdf.putVar(ncid, var_v_r, start, count, permute(v_r_slice, [2 1]));
    netcdf.putVar(ncid, var_v_i, start, count, permute(v_i_slice, [2 1]));
    
    delete(temp_file);
end

netcdf.close(ncid);
disp('Done. File created at:');
disp(output_file);
