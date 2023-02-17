#!/bin/bash
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -p gpu --gres=gpu:1
#SBATCH --mem=40G
#SBATCH --time 10:00:00
#SBATCH -J dwi_prepro
#SBATCH --output /gpfs/scratch/%u/dwi_prepro-log%A_%a.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hannah_swearingen1@brown.edu
#SBATCH --array=1-4


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

    ### Prepare files needed for topup ###
    ### Topup is not shell-dependent so it can be run in the subject loop to speed up processing ###
 

    ## PATH TO FIELDMAP/B0 FILES ##
    fieldmapfile_AP=${rawdata}/fmap/sub-${s}_ses-${ses}_acq-diffSE_dir-ap_epi.nii.gz
    fieldmapfile_PA=${rawdata}/fmap/sub-${s}_ses-${ses}_acq-diffSE_dir-pa_epi.nii.gz

    ## Create acquisition parameter and other files ##

    # Merge fieldmap/b0 files using < fslmerge -t outputfilename fieldmapfile_AP fieldmapfile_PA >

    fslmerge -t ${suboutdir}/AP_PA_b0 $fieldmapfile_AP $fieldmapfile_PA

    ## create acq_params text file for AP_PA_b0 file
    ## This acq_params file is used during topup
    ## This file should have a line for each volume in AP_PA_b0 but with respect to phase encoding direction and readout time.
    ## the acq_params text file is the same for each participant - so we will place it in the outdir

    if [ ! -f  "${outdir}/acq_params.txt" ]; then
	echo "......................................."
	echo "Creating acp_params text file for topup"
	echo "......................................."
	touch ${outdir}/acq_params.txt

	# number of volumes in fieldmap_AP
	ap_vols=$(fslinfo $fieldmapfile_AP | grep -m 1 dim4 | awk '{print $2}')
	pa_vols=$(fslinfo $fieldmapfile_PA | grep -m 1 dim4 | awk '{print $2}')

	for ((i=0; i<$ap_vols; i++)); do echo 0 1 0 ${readout_AP} >> ${outdir}/acq_params.txt; done
	for ((i=0; i<$pa_vols; i++)); do echo 0 -1 0 ${readout_PA} >> ${outdir}/acq_params.txt; done

    else
	echo ".............................................."
	echo "acq_params text file for topup already exists!"
	echo ".............................................."
    fi

    ########################### TOPUP #####################################

    ## Change into subject's directory in output dir ##

    cd $suboutdir
	  
    # Run topup. 
    # Output is AP_PA_topup_fieldcoef.nii.gz & AP_PA_topup_movpar.txt - files are only for subject/session not for each diffusion shell
    # fieldcoef are the inhomogeneity estimations
    # create additional output for brain mask

    echo "................................................................"
    echo "Running topup"
    echo "................................................................"

    topup --imain=AP_PA_b0.nii.gz --datain=${outdir}/acq_params.txt --config=b02b0.cnf --out=sub-${s}_ses-${ses}_AP_PA_topup --fout=sub-${s}_ses-${ses}_AP_PA_topup_HZ --iout=sub-${s}_ses-${ses}_topup_image

    echo "................................................................"
    echo "topup has completed for sub-$s"
    echo "................................................................"

    # Create a brain mask based on the first image from each topup -iout outup

    echo "................................................................"
    echo "Creating brain mask for sub-$s"
    echo "................................................................"

    fslroi sub-${s}_ses-${ses}_topup_image.nii.gz sub-${s}_ses-${ses}_firstvol.nii.gz 0 1
    bet sub-${s}_ses-${ses}_firstvol.nii.gz sub-${s}_ses-${ses}_brain.nii.gz -m -f 0.2

    echo "..............................................................."
    echo "Beginning loop through shells"
    echo "................................................................"    
   
    ### - Start shell loop - ###
    for shell in "${shells[@]}" 
    do
	
	echo "............................................................"
	echo "Defining variables for sub-${s} and ${shell} needed for eddy"
	echo "............................................................"

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

	######## EDDY ############

	if [ ! -f  "${outdir}/${shell}_acq_params.txt" ]; then
	    echo "..............................................."
	    echo "Creating ${shell}_acq_params text file for eddy"
	    echo "..............................................."

	    ## Create acq_params file - only needs to be two lines for eddy ##
	    touch ${outdir}/${shell}_acq_params.txt
	    echo 0 1 0  ${readout_AP} > ${outdir}/${shell}_acq_params.txt
	    echo 0 -1 0 ${readout_PA} >> ${outdir}/${shell}_acq_params.txt
	else
	    echo "............................................"
	    echo "${shell}_acq_params text file already exists!"
	    echo "............................................"
	fi

	if [ ! -f "${outdir}/${shell}_index.txt" ]; then
	    echo "..............................................."
	    echo "Creating ${shell}_index text file for eddy"
	    echo "..............................................."

	    ## Number of volumes in each diffusion scan ##
	    volumes_AP=$(fslinfo $diff_AP | grep -m 1 dim4 | awk '{print $2}')
	    volumes_PA=$(fslinfo $diff_PA | grep -m 1 dim4 | awk '{print $2}')

	    # Write the index file of the diffusion - needs to be a ROW
	    for i in $(eval echo "{1..$volumes_AP}"); do echo "1" >> ${outdir}/${shell}_column.txt; done
	    for i in $(eval echo "{1..$volumes_PA}"); do echo "2" >> ${outdir}/${shell}_column.txt; done
	    echo $(cat ${outdir}/${shell}_column.txt) | sed '/^$/d' > ${outdir}/${shell}_index.txt

        else
	    echo "............................................"
	    echo "${shell}_index text file already exists!"
	    echo "............................................"
	fi

        echo "................................................................"
	echo "Merging AP and PA diffusion, bvec, and bvals for sub-${s} ${shell}"
	echo "................................................................"

	# Merge AP and PA diffusion files #
	fslmerge -t ${suboutdir}/sub-${s}_ses-${ses}_acq-${shell}_combined_dwi $diff_AP $diff_PA

	# Combine AP and PA bvecs and bvals files #
	paste $bvecs_AP $bvecs_PA  > ${suboutdir}/bvecs_acq-${shell}.bvec
	paste $bvals_AP $bvals_PA  > ${suboutdir}/bvals_acq-${shell}.bval
	
	# For eddy options and directions for formatting sliceinfo.txt see:
	# https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddy/UsersGuide#A--index

	echo "................................................................"
	echo "Running eddy for sub-${s} ${shell}"
	echo "................................................................"

	eddy_cuda10.2 --imain=${suboutdir}/sub-${s}_ses-${ses}_acq-${shell}_combined_dwi.nii.gz \
		      --mask=sub-${s}_ses-${ses}_brain_mask.nii.gz \
		      --index=${outdir}/${shell}_index.txt \
		      --acqp=${outdir}/${shell}_acq_params.txt \
		      --bvecs=${suboutdir}/bvecs_acq-${shell}.bvec \
		      --bvals=${suboutdir}/bvals_acq-${shell}.bval \
		      --fwhm=2,1,0,0,0 \
		      --topup=sub-${s}_ses-${ses}_AP_PA_topup \
		      --out=sub-${s}_ses-${ses}_acq-${shell}_eddycorrected \
		      --repol \
		      --mporder=6 \
		      --json=$diff_json \
		      --estimate_move_by_susceptibility \
		      --verbose

	echo "................................................................"
	echo "Eddy has completed for sub-${s} ${shell}"
	echo "................................................................"

	echo "................................................................"
	echo "Running eddy QA for sub-${s} ${shell}"
	echo "................................................................"

	eddy_quad sub-${s}_ses-${ses}_acq-${shell}_eddycorrected \
		  -idx ${outdir}/${shell}_index.txt \
		  -par ${outdir}/${shell}_acq_params.txt \
		  -m ${suboutdir}/sub-${s}_ses-${ses}_brain_mask.nii.gz \
	          -b ${suboutdir}/bvals_acq-${shell}.bval \
		  -f sub-${s}_ses-${ses}_AP_PA_topup_HZ \
		  -o ${suboutdir}/sub-${s}_ses-${ses}_acq-${shell}_eddy_qc \
		  -v

	echo "................................................................"
	echo "Processing complete for sub-${s} ${shell}"
	echo "................................................................"

    done

    echo "................................................................"
    echo "Diffusion pre-processing complete for sub-${s}"
    echo "................................................................"
done
