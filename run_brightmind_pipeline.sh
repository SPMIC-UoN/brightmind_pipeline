#!/bin/bash

export SelfDir=$(dirname "$(readlink -f "$0")")

Usage()
{
  echo "Usage: `basename $0` [OPTIONS]"
  echo " "
  echo " "
  echo "OPTIONS:"
  echo " "
  echo " "
  echo " --path <path>                      (Mandatory) BRC Pipeline subjects path"
  echo " --subject <subject id>             (Mandatory) Subject ID"
  echo " --landmarks-file <file>            (Mandatory) Text file containing the coordinates in mm of the nasion and preauricular points"
  echo " --background-thresh <value>        (Mandatory) Integer specifying the intensity threshold between head and background on subject T1w image"
  echo " --output-path <path>               (Optional)  Path of the folder where the compressed directory with StimGuide files is written (Default: 'StimGuide' folder inside subject directory)"
  echo " --fmri-folder-name <folder name>   (Optional)  Name of the fMRI folder from the BRC functional analysis pipeline. Default: rfMRI"
  echo " -h | --help                        help"
  echo " "
  echo " "
}

# Just give usage if no arguments specified
if [ $# -eq 0 ] ; then Usage; exit 0; fi

# Mandatory arguments
Path=
Sub_ID=
Landmarks_File=
Bkgr_Thresh=

# Default values
Output_Path=
fMRI_Folder_Name="rfMRI"

# parse arguments
while [ "$1" != "" ]
do
    case $1 in
        --path )	            shift
                                Path=$1
                                ;;

        --subject )				shift
				                Sub_ID=$1
                                ;;
								
		--landmarks-file )		shift
                                Landmarks_File=$1
                                ;;
								
		--background-thresh )   shift
				                Bkgr_Thresh=$1
                                ;;
								
		--output-path )	        shift
                                Output_Path=$1
                                ;;
								
		--fmri-folder-name )	shift
                                fMRI_Folder_Name=$1
                                ;;
															
        -h | --help )           Usage
                                exit
                                ;;

        * )                     echo "Unknown flag: $1"
								echo
								Usage
                                exit 1
    esac
    shift
done

### Sanity checking of arguments
if [ X$Path = X ] || [ X$Sub_ID = X ] || [ X$Landmarks_File = X ] || [ X$Bkgr_Thresh = X ] ; then
  echo "All of the compulsory arguments --path, --subject, --landmarks-file and --background-thresh MUST be used"
  Usage
  exit 1;
fi

regex_integer='^[1-9][0-9]*$'
if ! [[ $Bkgr_Thresh =~ $regex_integer ]] ; then
	echo "Argument --background-thresh MUST be a valid integer"
	Usage
	exit 1;
fi

Scripts_Dir=$SelfDir/scripts
Root_Proc_Dir=$Path/$Sub_ID

if [ X$Output_Path = X ] ; then
	Output_Path=$Root_Proc_Dir/StimGuide/
fi

### Clean previous processed data (if it exists)
rm -rf $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg
rm -rf $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/granger
mkdir -p $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/melodic
mkdir -p $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/granger

echo -n "Preparing functional mask ... "
fslmaths $Root_Proc_Dir/analysis/anatMRI/T1/processed/data/T1_brain_mask.nii.gz -kernel sphere 3 -ero $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/func_mask.nii.gz -odt char
echo "done"

echo -n "Masking CSF and WM partial volume estimation images ... "
fslmaths $Root_Proc_Dir/analysis/anatMRI/T1/processed/seg/tissue/sing_chan/T1_pve_CSF.nii.gz -mas $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/func_mask.nii.gz $Root_Proc_Dir/analysis/anatMRI/T1/processed/seg/tissue/sing_chan/T1_pve_CSF_masked.nii.gz
fslmaths $Root_Proc_Dir/analysis/anatMRI/T1/processed/seg/tissue/sing_chan/T1_pve_WM.nii.gz -mas $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/func_mask.nii.gz $Root_Proc_Dir/analysis/anatMRI/T1/processed/seg/tissue/sing_chan/T1_pve_WM_masked.nii.gz
echo "done"

echo -n "Thresholding CSF and WM partial volume estimation images ... "
fslmaths $Root_Proc_Dir/analysis/anatMRI/T1/processed/seg/tissue/sing_chan/T1_pve_CSF_masked.nii.gz -sub 0.98 -bin $Root_Proc_Dir/analysis/anatMRI/T1/processed/seg/tissue/sing_chan/T1_pve_CSF_thresh.nii.gz
fslmaths $Root_Proc_Dir/analysis/anatMRI/T1/processed/seg/tissue/sing_chan/T1_pve_WM_masked.nii.gz -sub 0.98 -bin $Root_Proc_Dir/analysis/anatMRI/T1/processed/seg/tissue/sing_chan/T1_pve_WM_thresh.nii.gz
echo "done"

echo -n "Eroding CSF and WM partial volume estimation images ... "
fslmaths $Root_Proc_Dir/analysis/anatMRI/T1/processed/seg/tissue/sing_chan/T1_pve_CSF_thresh.nii.gz -kernel sphere 2 -ero $Root_Proc_Dir/analysis/anatMRI/T1/processed/seg/tissue/sing_chan/T1_pve_CSF_eroded.nii.gz
fslmaths $Root_Proc_Dir/analysis/anatMRI/T1/processed/seg/tissue/sing_chan/T1_pve_WM_thresh.nii.gz -kernel sphere 2 -ero $Root_Proc_Dir/analysis/anatMRI/T1/processed/seg/tissue/sing_chan/T1_pve_WM_eroded.nii.gz
echo "done"

echo -n "Transforming CSF and WM partial volume estimation images ... "
applywarp -i $Root_Proc_Dir/analysis/anatMRI/T1/processed/seg/tissue/sing_chan/T1_pve_CSF_eroded.nii.gz -r $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/data/SBref.nii.gz -o $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/CSF_highres2func.nii.gz -w $Root_Proc_Dir/analysis/$fMRI_Folder_Name/preproc/reg/str2rfMRI.nii.gz --interp=nn
applywarp -i $Root_Proc_Dir/analysis/anatMRI/T1/processed/seg/tissue/sing_chan/T1_pve_WM_eroded.nii.gz -r $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/data/SBref.nii.gz -o $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/WM_highres2func.nii.gz -w $Root_Proc_Dir/analysis/$fMRI_Folder_Name/preproc/reg/str2rfMRI.nii.gz --interp=nn
echo "done"

echo -n "Extracting WM and CSF timeseries ... "
fslmeants -i $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/data/rfMRI.nii.gz -o $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/WM_TS.txt -m $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/WM_highres2func.nii.gz
fslmeants -i $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/data/rfMRI.nii.gz -o $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/CSF_TS.txt -m $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/CSF_highres2func.nii.gz
echo "done"

echo -n "Creating regressor file ... "
paste $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/WM_TS.txt $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/CSF_TS.txt > $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/WM_CSF_TS.txt
echo "done"

echo -n "Regressing out WM and CSF from fMRI data ... "
fsl_regfilt -i $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/data/rfMRI.nii.gz -o $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/data/rfMRI_CSFWMregressed.nii.gz -d $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/WM_CSF_TS.txt -f "1,2"
echo "done"

echo -n "Transforming MNI F3 location into structural space ... "
std2imgcoord -img $Root_Proc_Dir/analysis/anatMRI/T1/processed/data/T1.nii.gz -std $SelfDir/MNI/MNI152_T1_1mm.nii.gz -xfm $Root_Proc_Dir/analysis/anatMRI/T1/preproc/reg/T1_2_std.mat -mm $SelfDir/MNI/F3_location.txt | awk '{printf("%.1f %.1f %.1f\n", $1, $2, $3)}' > $Root_Proc_Dir/analysis/anatMRI/T1/processed/data/F3_location.txt
echo "done"

echo -n "Transforming MNI masks into functional space ... "
applywarp -i $SelfDir/MNI/lCEN_front.nii.gz -r $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/data/SBref.nii.gz -o $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/lCEN_front_map.nii.gz -w $Root_Proc_Dir/analysis/$fMRI_Folder_Name/preproc/reg/std2rfMRI.nii.gz --interp=nn
applywarp -i $SelfDir/MNI/lMFG.nii.gz -r $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/data/SBref.nii.gz -o $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/lMFG_mask.nii.gz -w $Root_Proc_Dir/analysis/$fMRI_Folder_Name/preproc/reg/std2rfMRI.nii.gz --interp=nn
applywarp -i $SelfDir/MNI/rAI.nii.gz -r $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/data/SBref.nii.gz -o $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/rAI_mask.nii.gz -w $Root_Proc_Dir/analysis/$fMRI_Folder_Name/preproc/reg/std2rfMRI.nii.gz --interp=nn
echo "done"

echo -n "Running MELODIC to find independent components ... "
melodic -i $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/data/rfMRI_CSFWMregressed.nii.gz -o $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/melodic
echo "done"

echo -n "Thresholding and masking MELODIC components ... "
fslmaths $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/melodic/melodic_IC.nii.gz -thr 1.96 -mas $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/lCEN_front_map.nii.gz $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/melodic/melodic_IC_thresh.nii.gz
echo "done"

echo -n "Computing MELODIC component most correlated with Left Central Executive Network mask ... "
best_component=`fslcc -t 0.05 $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/lCEN_front_map.nii.gz $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/melodic/melodic_IC_thresh.nii.gz | sort -k3gr | head -1 | awk '{print $2}'`
let best_component--
echo "done"

echo -n "Extracting, binarising and masking MELODIC component (Number $best_component) ... "
fslroi $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/melodic/melodic_IC_thresh.nii.gz $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/melodic/melodic_IC_thresh_best.nii.gz $best_component 1
fslmaths $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/melodic/melodic_IC_thresh_best.nii.gz -mas $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/lMFG_mask.nii.gz -bin $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/melodic/melodic_IC_mask.nii.gz
echo "done"

echo -n "Running Granger analysis ... "
${MATLABpath}/matlab -nodesktop -nosplash -r "addpath(genpath('$Scripts_Dir/MATLAB')); compute_coefficient_GCA_maps('$Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/data/rfMRI_CSFWMregressed.nii.gz', '$Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/rAI_mask.nii.gz', '$Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/melodic/melodic_IC_mask.nii.gz', '$Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/granger/', 1); exit" > /dev/null
echo "done"

echo -n "Warping Granger map into structural space ... "
applywarp -i $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/granger/ZGCA_map_X2Y_1.nii.gz -r $Root_Proc_Dir/analysis/anatMRI/T1/processed/data/T1.nii.gz -o $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/granger/ZGCA_map_X2Y_1_T1.nii.gz -w $Root_Proc_Dir/analysis/$fMRI_Folder_Name/preproc/reg/rfMRI2str.nii.gz --interp=trilinear
echo "done"

echo -n "Warping MELODIC mask into structural space ... "
applywarp -i $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/melodic/melodic_IC_mask.nii.gz -r $Root_Proc_Dir/analysis/anatMRI/T1/processed/data/T1.nii.gz -o $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/melodic/melodic_IC_mask_T1.nii.gz -w $Root_Proc_Dir/analysis/$fMRI_Folder_Name/preproc/reg/rfMRI2str.nii.gz --interp=nn
echo "done"

echo -n "Computing Granger map statistics ... "
stats_out=`smoothest -z $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/granger/ZGCA_map_X2Y_1_T1.nii.gz -m $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/melodic/melodic_IC_mask_T1.nii.gz`
dlh=`echo $stats_out | awk '{print $2}'`
volume=`echo $stats_out | awk '{print $4}'`
echo "done"

# Find Granger target point (peak of most significant cluster at 95%. If not found, try at 90%, 80%, 70%, etc)
echo -n "Finding Granger target point ... "
cluster_peak=`cluster -i $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/granger/ZGCA_map_X2Y_1_T1.nii.gz -t 1.96 -p 0.05 -d ${dlh} --volume=${volume} --mm | tail -n +2 | sort -k3g | head -1 | awk '{printf("%.1f %.1f %.1f\n", $6, $7, $8)}'`
p_val="0.1"
while [[ "$cluster_peak" == "" ]]; do
	cluster_peak=`cluster -i $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/granger/ZGCA_map_X2Y_1_T1.nii.gz -t 1.96 -p ${p_val} -d ${dlh} --volume=${volume} --mm | tail -n +2 | sort -k3g | head -1 | awk '{printf("%.1f %.1f %.1f\n", $6, $7, $8)}'`
	p_val=`echo $p_val + 0.1 | bc`
done
echo $cluster_peak > $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/granger/${Sub_ID}_target_point_T1.txt
echo "done"

echo -n "Creating full landmark files ... "
cat -s $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/granger/${Sub_ID}_target_point_T1.txt $Landmarks_File > $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/data/landmarks_granger.txt
cat -s $Root_Proc_Dir/analysis/anatMRI/T1/processed/data/F3_location.txt $Landmarks_File > $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/data/landmarks_F3.txt
echo "done"

echo -n "Warping MNI head map and strip into structural space ... "
flirt -in $SelfDir/MNI/headmask_dil.nii.gz -ref $Root_Proc_Dir/analysis/anatMRI/T1/processed/data/T1.nii.gz -out $Root_Proc_Dir/analysis/anatMRI/T1/processed/data/T1_headmask_dil.nii.gz -applyxfm -init $Root_Proc_Dir/analysis/anatMRI/T1/preproc/reg/std_2_T1.mat -interp nearestneighbour
flirt -in $SelfDir/MNI/headstrip.nii.gz -ref $Root_Proc_Dir/analysis/anatMRI/T1/processed/data/T1.nii.gz -out $Root_Proc_Dir/analysis/anatMRI/T1/processed/data/T1_headstrip.nii.gz -applyxfm -init $Root_Proc_Dir/analysis/anatMRI/T1/preproc/reg/std_2_T1.mat -interp nearestneighbour
echo "done"

echo -n "Removing background of T1 image (threshold = $Bkgr_Thresh) ... "
fslmaths $Root_Proc_Dir/analysis/anatMRI/T1/processed/data/T1.nii.gz -mas $Root_Proc_Dir/analysis/anatMRI/T1/processed/data/T1_headmask_dil.nii.gz -thr $Bkgr_Thresh $Root_Proc_Dir/analysis/anatMRI/T1/processed/data/T1_nobackground.nii.gz
echo "done"

echo -n "Creating head surface mesh ... "
python2.7 $Scripts_Dir/Python/generate_surface_mesh.py $Sub_ID $Root_Proc_Dir/analysis/anatMRI/T1/processed/data/T1_nobackground.nii.gz --subjects-dir $Path --no-bias-corr > /dev/null
echo "done"

echo -n "Creating Granger StimGuide files ... "
python2.7 $Scripts_Dir/Python/generate_treatment_files.py $Sub_ID $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/data/landmarks_granger.txt $Output_Path/tmp_granger --subjects-dir $Path > $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/granger/target_info.txt
echo "done"

echo -n "Creating F3 StimGuide files ... "
python2.7 $Scripts_Dir/Python/generate_treatment_files.py $Sub_ID $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/data/landmarks_F3.txt $Output_Path/tmp_F3 --subjects-dir $Path > $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/granger/target_info.txt
echo "done"

echo -n "Creating compressed files ... "
zip -q -j $Output_Path/${Sub_ID}_cgiTBS.zip $Output_Path/tmp_granger/*
zip -q -j $Output_Path/${Sub_ID}_rTMS.zip $Output_Path/tmp_F3/*
rm -rf $Output_Path/tmp_granger
rm -rf $Output_Path/tmp_F3
echo "done"