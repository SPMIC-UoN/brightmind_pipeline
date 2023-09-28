# BRIGhTMIND Trial Processing Pipeline

## Table of Contents

* [Prerequisities](#prerequisities)
* [Script Usage](#script-usage)
* [Landmark Definition](#landmark-definition)
* [Background Threshold Definition](#background-threshold-definition)


<a id="prerequisities"></a>
## Prerequisities

The BRIGhTMIND Trial Processing Pipeline has the following software requirements:

1. A 64-bit Linux Operating System

2. [BRC Pipeline](https://github.com/SPMIC-UoN/BRC_Pipeline). During the trial, version [1.2.6](https://github.com/SPMIC-UoN/BRC_Pipeline/releases/tag/v1.2.6) was used and therefore is the version that should be used to replicate results.
	* The BRC Pipeline has also its own software requirements.

3. MATLAB R2017b or later.

4. Python version 2.7 (the scripts should in principle also work with modern versions of Python (e.g. Python 3.x), but this has not been tested).

-----

<a id="script-usage"></a>
## Script Usage

Before running this script, the structural and resting-state fMRI datasets of the subject to process must be run through the BRC Pipeline's structural `BRC_structural_pipeline/struc_preproc.sh` and functional `BRC_functional_pipeline/fMRI_preproc.sh` procedures, respectively (in that order).

The following calls to BRC Pipeline were used in the BRIGhTMIND study:

* Structural pipeline

		struc_preproc.sh --input <path_to_T1w_nifti> \
		                 --path <brc_output_path> \
		                 --subject <subject id> \
		                 --strongbias \
		                 --nodefacing
	
* Functional pipeline

		fMRI_preproc.sh --fmripath <path_to_rsfMRI_nifti> \
		                --path <brc_output_path> \
		                --subject <subject id> \
		                --mctype MCFLIRT6 \
		                --dcmethod TOPUP \
		                --slice2vol \
		                --slspec <dcm2niix_json_file_for_rsfMRI> \
		                --fwhm 5 \
		                --echospacing <effective_echo_spacing_of_acquisition> \
		                --stcmethod 1 \
		                --SEPhasePos <path_to_rsfMRI_anterior_to_posterior_blip_nifti> \
		                --SEPhaseNeg <path_to_rsfMRI_posterior_to_anterior_blip_nifti> \
		                --unwarpdir y- \
		                --biascorrection SEBASED \
		                --echodiff NONE \
		                --fmapmag NONE \
		                --tempfilter 100
						
-----

The usage of run_brightmind_pipeline.sh is as follows:

         run_brightmind_pipeline.sh <arguments>
		 
* Compulsory arguments (You MUST set all of):

        --path <brc_output_path>           BRC Pipeline output path
        --subject <subject id>             Subject ID as set in BRC Pipeline
		--landmarks-file <file path>       Text file containing the coordinates in mm of the nasion and preauricular points
		--background-thresh <value>        Integer specifying the intensity threshold between head and background on subject T1w image

* Optional arguments (You may optionally specify one or more of):

        --output-path <path>               Path of the folder where the compressed directories with StimGuide files are written (Default: 'StimGuide' folder inside subject directory)
        --fmri-folder-name <folder name>   Name of the fMRI folder from the BRC functional analysis pipeline. Default: 'rfMRI'"
        -h | --help                        help
		
-----

<a id="landmark-definition"></a>
## Landmark Definition

The BRIGhTMIND processing pipeline requires as input a text file with the location in mm of the nasion, left preauricular point and right preauricular point. In order to create this file, open the file `<brc_output_path>/<subject id>/analysis/anatMRI/T1/processed/data/T1.nii.gz` in your favourite NIFTI image viewer.

### Location of the nasion

Locate the nasion as shown in the image below and take note of the X, Y and Z coordinates **in mm**

![Nasion Landmark](https://github.com/SPMIC-UoN/brightmind_pipeline/blob/main/img/nasion_landmark.png)

### Location of the left preauricular point

Locate the left preauricular point as shown in the image below and take note of the X, Y and Z coordinates **in mm**

![Left Preauricular Landmark](https://github.com/SPMIC-UoN/brightmind_pipeline/blob/main/img/left_preauricular_landmark.png)

### Location of the right preauricular point

Locate the right preauricular point as shown in the image below and take note of the X, Y and Z coordinates **in mm**

![Right Preauricular Landmark](https://github.com/SPMIC-UoN/brightmind_pipeline/blob/main/img/right_preauricular_landmark.png)

### Landmark File

The final landmark file should look like the image below. The first line corresponds to the nasion coordinates, the second line to the left preauricular point coordinates and the third line to the right preauricular point coordinates.

![Landmark Coordinates](https://github.com/SPMIC-UoN/brightmind_pipeline/blob/main/img/landmark_coords.png)

-----

<a id="background-threshold-definition"></a>
## Background Threshold Definition

Another mandatory input to the BRIGhTMIND pipeline is an integer number defining the threshold between background and foreground in the T1-weighted structural image. In general terms, the procedure to determine this number is as follows:

1. Open the file `<brc_output_path>/<subject id>/analysis/anatMRI/T1/processed/data/T1.nii.gz` in your favourite NIFTI image viewer.

2. Window the image by lowering the maximum intensity, such that the skull/scalp has a clear boundary with the background, but at the same time there is not too much noise near the skull/scalp like in the example below:

	![Background Threshold Setting](https://github.com/SPMIC-UoN/brightmind_pipeline/blob/main/img/background_thresh_setting.png)
	
3. The maximum intensity selected is the number to input into the script.

After the pipeline is run, the resulting GIFTI mesh `<brc_output_path>/<subject id>/surf/lh.seghead.surf.gii` should look smooth and complete like in the example below:

![Smooth Scalp Mesh](https://github.com/SPMIC-UoN/brightmind_pipeline/blob/main/img/scalp_smooth.png)

If the threshold number is set too low, the resulting scalp mesh may become too noisy like in this example:

![Rough Scalp Mesh](https://github.com/SPMIC-UoN/brightmind_pipeline/blob/main/img/scalp_rough.png)

If the threshold number is set too high, the resulting scalp mesh may be not complete and have holes in it.
