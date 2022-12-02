#!/bin/bash
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -p gpu --gres=gpu:1
#SBATCH --mem=40G
#SBATCH --time 10:00:00
#SBATCH -J dwi_prepro
#SBATCH --output /gpfs/scratch/%u/dwi_prepro-log%A_%a.txt
#SBATCH --array=1


###### CONFIGURE MODULES #######

module load fsl/6.0.5.2
module load cuda/10.2 cudnn/8.2.0

# This script is built on FSL and is adapted from the HCP-A processing scripts and assumes BIDS organization.
# It can be run on dwi collected on a sphere and requires a merged file of b0 images for eddy correction.
# This preprocessing approach is not appropriate for grid dwi acquisitions.

### - IMPORTANT MANUAL USER INPUT - ###
# The topup function needs an acquisition parameters (acq_params.txt) text file describing your spin echo field maps or b0 image (e.g. bo_AP and b0_PA or acq-diffSE_dir-ap/acq-diffSE_dir-pa).
# Examples can be found at: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/topup/TopupUsersGuide#A--datain
# The first three numbers in each row describe how image was phase encoded (x-, y-, z-dimensions; 'blip up, blip down')).
# The last number desribe the read-out time, in milliseconds.
# Extract this from the “TotalReadOutTime” field of the b0 or fieldmap .json file.
# You should have one row for each volume in your merged file.

readout_AP=0.104248
readout_PA=0.104248


### - Define variables and paths - ###

ses=01
codedir=/gpfs/data/jbarredo/diffusion/code #path to code directory
rawdir=/gpfs/data/jbarredo/diffusion/rawdata #path to rawdata directory
derivdir=/gpfs/data/jbarredo/diffusion/derivatives #path to derivatives directory
outdir=$derivdir/dwi_prepro  #path to dmri directory
subjects="`sed -n ${SLURM_ARRAY_TASK_ID}p dmri_subjects.txt`"
#subfile=dmri_subjects.txt #text file with list of subjects
shells=(b1500) #array for shells (1000 2000)

### IMPORTANT - IF YOU CHANGE FROM USING DIRECTIONS TO RUNS OR REMOVE THE DIRECTIONS VARIABLE YOU WILL NEED TO UPDATE THE SCRIPT TO REFLECT YOUR DATA'S NAMING CONVENTION OR REMOVE THE DIRECTIONS FOR LOOP ###

### - Start subject loop - ###
### Subjects loop will loop through each subject in subjects file and create directory variables ###

for s in $subjects #loop through subjects array
do
    rawdata=${rawdir}/sub-${s}/ses-${ses}
    mkdir ${outdir}/sub-${s}
    mkdir ${outdir}/sub-${s}/ses-${ses}
    suboutdir=${outdir}/sub-${s}/ses-${ses}

    ### - Start shell loop - ###
    ### Shell loop will combine AP/PA diffusion files, bvec & bval files, and run applytopup and eddy after the directions loop runs topup ###

    for shell in "${shells[@]}" 
    do
	## PATH TO DIFFUSION FILES ##
	diff_AP=${rawdata}/dwi/sub-${s}_ses-${ses}_acq-${shell}_dir-ap_dwi.nii.gz
	diff_PA=${rawdata}/dwi/sub-${s}_ses-${ses}_acq-${shell}_dir-pa_dwi.nii.gz

	# path to a json file - eddy will extract the slice times from the json file. Only one json file is needed so either ap/pa will work #
	diff_json=${rawdata}/dwi/sub-${s}_ses-${ses}_acq-${shell}_dir-ap_dwi.json

	## Set bvecs and bvals variables ##
	bvecs_AP=${rawdata}/dwi/sub-${s}_ses-${ses}_acq-${shell}_dir-ap_dwi.bvec
	bvals_AP=${rawdata}/dwi/sub-${s}_ses-${ses}_acq-${shell}_dir-ap_dwi.bval

      	bvecs_PA=${rawdata}/dwi/sub-${s}_ses-${ses}_acq-${shell}_dir-pa_dwi.bvec
	bvals_PA=${rawdata}/dwi/sub-${s}_ses-${ses}_acq-${shell}_dir-pa_dwi.bval

	## PATH TO FIELDMAP/B0 FILES ##
	fieldmapfile_AP=${rawdata}/fmap/sub-${s}_ses-${ses}_acq-diffSE_dir-ap_epi.nii.gz
	fieldmapfile_PA=${rawdata}/fmap/sub-${s}_ses-${ses}_acq-diffSE_dir-pa_epi.nii.gz

	## Create acquisition parameter and other files ##

	# Merge fieldmap/b0 files using < fslmerge -t outputfilename fieldmapfile_AP fieldmapfile_PA >

	fslmerge -t ${suboutdir}/AP_PA_b0 $fieldmapfile_AP $fieldmapfile_PA

	## create acq_params text file for AP_PA_b0 file
	## This acq_params file is used during topup
	## This file should have a line for each volume in AP_PA_b0 but with respect to phase encoding direction and readout time.

	touch ${suboutdir}/acq_params.txt

	# number of volumes in fieldmap_AP
	ap_vols=$(fslinfo $fieldmapfile_AP | grep -m 1 dim4 | awk '{print $2}')
	pa_vols=$(fslinfo $fieldmapfile_PA | grep -m 1 dim4 | awk '{print $2}')

	for ((i=0; i<$ap_vols; i++)); do echo 0 1 0 ${readout_AP} >> ${suboutdir}/acq_params.txt; done
	for ((i=0; i<$pa_vols; i++)); do echo 0 -1 0 ${readout_PA} >> ${suboutdir}/acq_params.txt; done
	    
	#the acq_params.txt will be removed at each iteration or else it will continually be appended  

	  
	########################### TOPUP #####################################

	 ## Change into subject's directory in output dir ##

	cd $suboutdir
	  
	# Run topup. 
	# Output is AP_PA_topup_fieldcoef.nii.gz & AP_PA_topup_movpar.txt - files are only for subject/session not for each diffusion shell
	# fieldcoef are the inhomogeneity estimations

	topup --imain=AP_PA_b0.nii.gz --datain=${suboutdir}/acq_params.txt --config=b02b0.cnf --out=sub-${s}_ses-${ses}_AP_PA_topup --fout=sub-${s}_ses-${ses}_AP_PA_topup_HZ


	# Apply correction to the diffusion files (separate each diffusion direction by comma for --imain)
	# The index should correspond to lines of acq_param.txt from the same phase-encoding direction
	# Output will be a topup_corrected file corresponding to each shell (if applicable)

	applytopup --imain=$diff_AP,$diff_PA \
	    --inindex=1,2 --datain=${suboutdir}/acq_params.txt \
	    --topup=${suboutdir}/sub-${s}_ses-${ses}_AP_PA_topup --method=jac \
	    --verbose --out=${suboutdir}/sub-${s}_ses-${ses}_${shell}_topup_corrected

	######## EDDY ############

	## Number of volumes in each diffusion scan ##
	volumes_AP=$(fslinfo $diff_AP | grep -m 1 dim4 | awk '{print $2}')
	volumes_PA=$(fslinfo $diff_PA | grep -m 1 dim4 | awk '{print $2}')

	# Write the index file of the diffusion - needs to be a ROW
	for i in $(eval echo "{1..$volumes_AP}"); do echo "1" >> ${suboutdir}/${shell}_column.txt; done #removed after each iteration
	for i in $(eval echo "{1..$volumes_PA}"); do echo "2" >> ${suboutdir}/${shell}_column.txt; done
	echo $(cat ${suboutdir}/${shell}_column.txt) | sed '/^$/d' > ${suboutdir}/${shell}_index.txt #removed after each iteration

	## Create acq_params file - only needs to be two lines for eddy ##
	touch ${suboutdir}/${shell}_acq_params.txt
	echo 0 1 0  ${readout_AP} > ${suboutdir}/${shell}_acq_params.txt
	echo 0 -1 0 ${readout_PA} >> ${suboutdir}/${shell}_acq_params.txt

	# Create a brain mask based on the first image from each topup corrected dataset.
	fslroi sub-${s}_ses-${ses}_${shell}_topup_corrected.nii.gz sub-${s}_ses-${ses}_${shell}_firstvol_corrected.nii.gz 0 1
	bet sub-${s}_ses-${ses}_${shell}_firstvol_corrected.nii.gz sub-${s}_ses-${ses}_${shell}_brain.nii.gz -m -f 0.2

	# Merge AP and PA diffusion files #
	fslmerge -t sub-${s}_ses-${ses}_${shell}_combined_dwi $diff_AP $diff_PA

	# Combine AP and PA bvecs and bvals files #
	paste $bvecs_AP $bvecs_PA  > bvecs_all.bvec
	paste $bvals_AP $bvals_PA  > bvals_all.bval
	
	# For eddy options and directions for formatting sliceinfo.txt see:
	# https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddy/UsersGuide#A--index

	eddy_cuda10.2 --imain=sub-${s}_ses-${ses}_${shell}_combined_dwi.nii.gz \
		      --mask=sub-${s}_ses-${ses}_${shell}_brain_mask.nii.gz \
		      --index=${suboutdir}/${shell}_index.txt \
		      --acqp=${suboutdir}/${shell}_acq_params.txt \
		      --bvecs=${suboutdir}/bvecs_all.bvec \
		      --bvals=${suboutdir}/bvals_all.bval \
		      --fwhm=2,1,0,0,0 \
		      --topup=sub-${s}_ses-${ses}_AP_PA_topup \
		      --out=sub-${s}_ses-${ses}_${shell}_eddycorrected \
		      --repol \
		      --mporder=6 \
		      --json=$diff_json \
		      --estimate_move_by_susceptibility \
		      --verbose

	eddy_quad sub-${s}_ses-${ses}_${shell}_eddycorrected \
		  -idx ${suboutdir}/${shell}_index.txt \
		  -par ${suboutdir}/${shell}_acq_params.txt \
		  -m sub-${s}_ses-${ses}_${shell}_brain_mask.nii.gz \
	          -b ${suboutdir}/bvals_all.bval \
		  -f sub-${s}_ses-${ses}_AP_PA_topup_HZ \
		  -o ${suboutdir}/sub-${s}_ses-${ses}_${shell}_eddy_qc \
		  -v
    done

    rm ${suboutdir}/${shell}_index.txt
    rm ${suboutdir}/${shell}_column.txt
    rm ${suboutdir}/${shell}_acq_params.txt
    rm ${suboutdir}/acq_params.txt

done
