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
#SBATCH --array=1-57

# This script is built on FSL and is adapted from the HCP-A processing scripts and assumes BIDS organization.
# It can be run on dwi collected on a sphere and requires a merged file of b0 images for eddy correction.
# This preprocessing approach is not appropriate for grid dwi acquisitions.
# This script is for dwi acquisitions NOT collected in different phase-encoding directions

###### CONFIGURE MODULES #######

module load fsl/6.0.5.2
module load cuda/10.2 cudnn/8.2.0

### - IMPORTANT MANUAL USER INPUT - ###
# The topup function needs an acquisition parameters (acq_params.txt) text file describing your spin echo field maps or b0 image (e.g. bo_AP and b0_PA or acq-diffSE_dir-ap/acq-diffSE_dir-pa).
# Examples can be found at: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/topup/TopupUsersGuide#A--datain
# The first three numbers in each row describe how image was phase encoded (x-, y-, z-dimensions; 'blip up, blip down')).
# The last number desribe the read-out time, in milliseconds.
# Extract this from the “TotalReadOutTime” field of the b0 or fieldmap .json file.
# You should have one row for each volume in your merged file.


readout_AP=0.0334748
readout_PA=0.0345046


### - Define variables and paths - ###

ses=01
codedir=/gpfs/data/jbarredo/nstb/code #path to code directory
rawdir=/gpfs/data/jbarredo/nstb/rawdata #path to rawdata directory
derivdir=/gpfs/data/jbarredo/nstb/derivatives #path to derivatives directory
outdir=$derivdir/dmriprepro  #path to dmri directory
subjects="`sed -n ${SLURM_ARRAY_TASK_ID}p dmri_subjects.txt`"
#subfile=dmri_subjects.txt #text file with list of subjects - file created in code directory
shells=(1000 2000 3000) #array for shells (1000 2000)

### - Start subject loop - ###
### Subjects loop will loop through each subject in subjects file and create directory variables ###

for s in $subjects # loop through subjects array
do
    rawdata=${rawdir}/sub-${s}/ses-${ses}/dwi
    mkdir ${outdir}/sub-${s}
    mkdir ${outdir}/sub-${s}/ses-${ses}
    suboutdir=${outdir}/sub-${s}/ses-${ses}

    ### Prepare files needed for topup ###
    ### Topup is not shell-dependent so it can be run in the subject loop to speed up processing ###

    ## PATH TO FIELDMAP/B0 FILES ##
    fieldmapfile_AP=${rawdata}/sub-${s}_ses-${ses}_acq-b0_dir-AP_dwi.nii.gz
    fieldmapfile_PA=${rawdata}/sub-${s}_ses-${ses}_acq-b0_dir-PA_dwi.nii.gz

    ## Create merged b0 file, acquisition parameter file, and other text files ##

    # create the file “AP_PA_b0.nii.gz” from b0_AP and b0_PA
    fslmerge -t ${suboutdir}/AP_PA_b0 $fieldmapfile_AP $fieldmapfile_PA

    ##create B0_acq_params text file
    ## This acq_params file is used during topup
    ## This file should have a line for each volume in AP_PA_b0 but with respect to phase encoding direction and readout time.
    touch ${suboutdir}/acq_params.txt

    # number of volumes in fieldmap_AP
    ap_vols=$(fslinfo $fieldmapfile_AP | grep -m 1 dim4 | awk '{print $2}')
    pa_vols=$(fslinfo $fieldmapfile_PA | grep -m 1 dim4 | awk '{print $2}')

    for ((i=0; i<$ap_vols; i++)); do echo 0 1 0 ${readout_AP} >> ${suboutdir}/acq_params.txt; done
    for ((i=0; i<$pa_vols; i++)); do echo 0 -1 0 ${readout_PA} >> ${suboutdir}/acq_params.txt; done
	   
    ########################### TOPUP #####################################

    ## Change into subject's directory in output dir ##
    cd $suboutdir	
	
    # Run topup. Output is AP_PA_topup_fieldcoef.nii.gz & AP_PA_topup_movpar.txt
    # fieldcoef are the inhomogeneity estimations
    topup --imain=AP_PA_b0.nii.gz --datain=${suboutdir}/acq_params.txt --config=b02b0.cnf --out=sub-${s}_ses-${ses}_AP_PA_topup --fout=sub-${s}_ses-${ses}_AP_PA_topup_HZ  

    ### - Start shell loop - ###

    for shell in "${shells[@]}" # path to shells file
    do
	# PATH TO DIFFUSION FILES ##
	diff_file=${rawdata}/sub-${s}_ses-${ses}_acq-${shell}_dwi.nii.gz
	diff_json=${rawdata}/sub-${s}_ses-${ses}_acq-${shell}_dwi.json

	## Set additional variables ##
	bvecs=${rawdata}/sub-${s}_ses-${ses}_acq-${shell}_dwi.bvec
	bvals=${rawdata}/sub-${s}_ses-${ses}_acq-${shell}_dwi.bval	

	## Number of volumes in each diffusion scan ##
	volumes=$(fslinfo ${rawdata}/sub-${s}_ses-${ses}_acq-${shell}_dwi.nii.gz | grep -m 1 dim4 | awk '{print $2}')

	# Write the index file of the diffusion file - needs to be a ROW
	# Use "1" for AP direction; Use "2" for PA direction
	for i in $(eval echo "{1..$volumes}"); do echo "1" >> ${suboutdir}/${shell}_column.txt; done 
	echo $(cat ${suboutdir}/${shell}_column.txt) | sed '/^$/d' > ${suboutdir}/${shell}_index.txt

	## Create acq_params file for specific diffusion files ##
	## If diffusion file is only collected in 1 direction, only need 1 line for this file ##
	touch ${suboutdir}/${shell}_acq_params.txt
	echo 0 1 0  ${readout_AP} > ${suboutdir}/${shell}_acq_params.txt # AP encoding

	########## APPLYTOPUP ###########

	cd $suboutdir

	# Apply correction to the diffusion files. The index should correspond to lines of acq_param.txt from the same phase-encoding direction. 
	applytopup --imain=$diff_file \
		   --inindex=1 --datain=${suboutdir}/${shell}_acq_params.txt \
		   --topup=sub-${s}_ses-${ses}_AP_PA_topup --method=jac \
		   --verbose --out=${suboutdir}/sub-${s}_ses-${ses}_${shell}_topup_corrected

	######## EDDY ############

	# Create a brain mask based on the first image from each topup corrected dataset.
	fslroi sub-${s}_ses-${ses}_${shell}_topup_corrected.nii.gz sub-${s}_ses-${ses}_${shell}_firstvol_corrected.nii.gz 0 1
	bet sub-${s}_ses-${ses}_${shell}_firstvol_corrected.nii.gz sub-${s}_ses-${ses}_${shell}_brain.nii.gz -m -f 0.2

	# For eddy options and directions for formatting sliceinfo.txt see:
	# https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddy/UsersGuide#A--index

	eddy_cuda10.2 --imain=$diff_file \
		      --mask=sub-${s}_ses-${ses}_${shell}_brain_mask.nii.gz \
		      --index=${suboutdir}/${shell}_index.txt \
		      --acqp=${suboutdir}/${shell}_acq_params.txt \
		      --bvecs=$bvecs \
		      --bvals=$bvals \
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
		  -b $bvals \
		  -f sub-${s}_ses-${ses}_AP_PA_topup_HZ \
		  -o ${suboutdir}/sub-${s}_ses-${ses}_${shell}_eddy_qc \
		  -v
    done
done

