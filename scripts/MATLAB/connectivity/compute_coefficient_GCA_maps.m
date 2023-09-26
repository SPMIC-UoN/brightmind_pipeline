function compute_coefficient_GCA_maps(nii_fmri_in_path, nii_seed_mask_in_path, nii_searchlight_mask_in_path, output_dir, order)
%
%  compute_coefficient_GCA_maps(nii_fmri_in_path, nii_seed_mask_in_path, nii_searchlight_mask_in_path, output_dir, order)
%
%  Computes coefficient-based Granger Causality Analysis (GCA) maps for each voxel in a given searchlight with respect to a given seed.
%
%  Requires MATLAB R2017b or later.
%
%  Required MATLAB library:
%
%  Tools for NIfTI and ANALYZE image (https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image).
%
%  Input MATLAB variables:
%
%    - nii_fmri_in_path               : (Required) Path to the Nifti file of the fMRI timeseries, as a string or character array.
%    - nii_seed_mask_in_path          : (Required) Path to the Nifti file of the seed mask, as a string or character array. Voxels with a value
%                                       greater than zero are considered as being inside the mask.
%    - nii_searchlight_mask_in_path   : (Required) Path to the Nifti file of the searchlight mask, as a string or character array. Voxels with
%                                       a value greater than zero are considered as being inside the mask.
%    - output_dir                     : (Required) Directory of output Nifti files, as a string or character array.
%    - order                          : (Optional) Integer value specifying the time order of Granger analysis, i.e. how many timepoints (TR) in the 
%                                       past we look in the autogregressive models. Default is 1. 
%
%
%  Resulting Nifti files:
%
%    - GCA_map_X2Y_{P}.nii.gz  : Nifti file containing the seed to target Granger map for the time lag of P timepoints (P = 1 ... order).
%    - GCA_map_Y2X_{P}.nii.gz  : Nifti file containing the target to seed Granger map for the time lag of P timepoints (P = 1 ... order).
%    - ZGCA_map_X2Y_{P}.nii.gz : Nifti file containing the seed to target Z-scored Granger map for the time lag of P timepoints (P = 1 ... order).
%    - ZGCA_map_Y2X_{P}.nii.gz : Nifti file containing the target to seed Z-scored Granger map for the time lag of P timepoints (P = 1 ... order).
%
%                       
%  (c) 2023-2024. Written by Dr. Stefan Pszczolkowski (mszspp@nottingham.ac.uk) 
%
    if nargin < 4
        error('Not enough input arguments.');
    elseif nargin == 4
        order = 1;
    end
    
    output_dir = strrep(output_dir, '\', '/');
    prefix_path = fileparts(output_dir);
    if ~isempty(prefix_path) && ~exist(prefix_path, 'dir')
        error('Output folder ''%s'' does not exist.',prefix_path); 
    end
    
    if ~exist(nii_fmri_in_path, 'file')
        error('File ''%s'' does not exist.', nii_fmri_in_path); 
    end
    
    if ~exist(nii_seed_mask_in_path, 'file')
        error('File ''%s'' does not exist.', nii_seed_mask_in_path); 
    end
    
    if ~exist(nii_searchlight_mask_in_path, 'file')
        error('File ''%s'' does not exist.', nii_searchlight_mask_in_path); 
    end
    
    if ~isscalar(order) || ~isnumeric(order) || (rem(order, 1) ~= 0) || (order <= 0)
        error('Order must be a positive integer.'); 
    end
    
    % Make sure string-type input are (converted to) character arrays
    output_dir = convertStringsToChars(output_dir);
    nii_fmri_in_path = convertStringsToChars(nii_fmri_in_path);
    nii_seed_mask_in_path = convertStringsToChars(nii_seed_mask_in_path);
    nii_searchlight_mask_in_path = convertStringsToChars(nii_searchlight_mask_in_path);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    nii_fmri = load_untouch_nii(nii_fmri_in_path);
    nii_mask_seed = load_untouch_nii(nii_seed_mask_in_path);
    nii_mask_searchlight = load_untouch_nii(nii_searchlight_mask_in_path);
    
    if ndims(nii_fmri.img) ~= 4
        error('Provided fMRI image must be 4-dimensional');
    end
    
    [timeseries_size_x, timeseries_size_y, timeseries_size_z, timeseries_size_t] = size(nii_fmri.img);
    
    if ~isequal(size(nii_mask_seed.img), [timeseries_size_x, timeseries_size_y, timeseries_size_z])
       error('The size of the seed mask image does not match the volume size of the fMRI timeseries.'); 
    end
    
    if ~isequal(size(nii_mask_searchlight.img), [timeseries_size_x, timeseries_size_y, timeseries_size_z])
       error('The size of the searchlight mask image does not match the volume size of the fMRI timeseries.'); 
    end

    auto_regressive_data_length = timeseries_size_t - order;
    volume_size = [timeseries_size_x, timeseries_size_y, timeseries_size_z];
    volume_length = prod(volume_size);

    timepoint_data = reshape(nii_fmri.img, [volume_length, timeseries_size_t]);
    timepoint_data = normalize(timepoint_data, 2, 'center');
    
    score_map_x2y = zeros(volume_length, order);
    score_map_y2x = zeros(volume_length, order);
	score_map_x2y_zscores = zeros(volume_length, order);
    score_map_y2x_zscores = zeros(volume_length, order);
    
    mask_brain = (var(timepoint_data, [], 2) > 0);
    
    mask_seed = (nii_mask_seed.img(:) > 0) & mask_brain;
    mask_searchlight = (nii_mask_searchlight.img(:) > 0) & mask_brain & (~mask_seed);
    
    X = mean(timepoint_data(mask_seed,:))';
    X_reg = X(order+1:end);
    
    auto_regressive_data_x = zeros(auto_regressive_data_length, order);
    
    for t = 1:auto_regressive_data_length
        auto_regressive_data_x(t,:) = X(t:t+order-1)';
    end

    searchlight_indices = find(mask_searchlight)';
    
    for index = searchlight_indices   
        Y = timepoint_data(index,:)';
        Y_reg = Y(order+1:end);

        auto_regressive_data_y = zeros(auto_regressive_data_length, order);

        for t = 1:auto_regressive_data_length
            auto_regressive_data_y(t,:) = Y(t:t+order-1)';
        end

        auto_regressive_data_combined = [auto_regressive_data_x, auto_regressive_data_y];

        beta_x2y = regress(Y_reg, auto_regressive_data_combined);
        beta_y2x = regress(X_reg, auto_regressive_data_combined);
    
        score_map_x2y(index,:) = beta_x2y(1:order);
        score_map_y2x(index,:) = beta_y2x(order+1:2*order);
    end

	for p = 1:order
		score_map_x2y_zscores(:,p) = (score_map_x2y(:,p) - mean(score_map_x2y(:,p))) ./ std(score_map_x2y(:,p));
		score_map_y2x_zscores(:,p) = (score_map_y2x(:,p) - mean(score_map_y2x(:,p))) ./ std(score_map_y2x(:,p));
	end

    for p = 1:order
        nii_x2y_out = nii_mask_seed;
        nii_x2y_out.img = reshape(score_map_x2y(:,p), volume_size);
        nii_x2y_out.hdr.dime.bitpix = 32;
        nii_x2y_out.hdr.dime.datatype = 16;
        nii_x2y_out.hdr.dime.glmin = min(nii_x2y_out.img(:));
        nii_x2y_out.hdr.dime.glmax = max(nii_x2y_out.img(:));
        save_untouch_nii(nii_x2y_out, [output_dir '/GCA_map_X2Y_' num2str(p) '.nii.gz']);

        nii_y2x_out = nii_mask_seed;
        nii_y2x_out.img = reshape(score_map_y2x(:,p), volume_size);
        nii_y2x_out.hdr.dime.bitpix = 32;
        nii_y2x_out.hdr.dime.datatype = 16;
        nii_y2x_out.hdr.dime.glmin = min(nii_y2x_out.img(:));
        nii_y2x_out.hdr.dime.glmax = max(nii_y2x_out.img(:));
        save_untouch_nii(nii_y2x_out, [output_dir '/GCA_map_Y2X_' num2str(p) '.nii.gz']);
		
		nii_x2y_zscores_out = nii_mask_seed;
        nii_x2y_zscores_out.img = reshape(score_map_x2y_zscores(:,p), volume_size);
        nii_x2y_zscores_out.hdr.dime.bitpix = 32;
        nii_x2y_zscores_out.hdr.dime.datatype = 16;
        nii_x2y_zscores_out.hdr.dime.glmin = min(nii_x2y_zscores_out.img(:));
        nii_x2y_zscores_out.hdr.dime.glmax = max(nii_x2y_zscores_out.img(:));
        save_untouch_nii(nii_x2y_zscores_out, [output_dir '/ZGCA_map_X2Y_' num2str(p) '.nii.gz']);

        nii_y2x_zscores_out = nii_mask_seed;
        nii_y2x_zscores_out.img = reshape(score_map_y2x_zscores(:,p), volume_size);
        nii_y2x_zscores_out.hdr.dime.bitpix = 32;
        nii_y2x_zscores_out.hdr.dime.datatype = 16;
        nii_y2x_zscores_out.hdr.dime.glmin = min(nii_y2x_zscores_out.img(:));
        nii_y2x_zscores_out.hdr.dime.glmax = max(nii_y2x_zscores_out.img(:));
        save_untouch_nii(nii_y2x_zscores_out, [output_dir '/ZGCA_map_Y2X_' num2str(p) '.nii.gz']);
    end
end

