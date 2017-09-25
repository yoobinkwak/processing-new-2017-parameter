FSLDIR=/usr/share/fsl/5.0
FREESURFER_HOME=/usr/local/freesurfer
export SUBJECTS_DIR=/home/yoobinkwak/preproc_BCS_001_002
HCP_scripts=Pipelines-master/FreeSurfer/scripts

for i in $@
do
    FS=${i}/Preprocess_Structure/Freesurfer
    if [ ! -d ${FS} ]
    then
        mkdir ${FS}
    fi

    T1=${i}/Preprocess_Structure/PreFreesurfer/OutputT1wRestoredImage.nii.gz
    T1_brain=${i}/Preprocess_Structure/PreFreesurfer/OutputT1wRestoredBrainImage.nii.gz
    T2=${i}/Preprocess_Structure/PreFreesurfer/OutputT2wRestoredImage.nii.gz

#### downsample T1 (bias corrected) to 1mm ####
    downsampled=${FS}/1mm_downsampled
    if [ ! -d ${downsampled} ]
    then
        mkdir ${downsampled}
    fi
    T1_1mm=${downsampled}/T1_1mm.nii.gz
    T1_brain_1mm=${downsampled}/T1_brain_1mm.nii.gz
    if [ ! -f ${T1_1mm} ]
    then
        Mean=`${FSLDIR}/bin/fslstats ${T1_brain} -M`
        ${FSLDIR}/bin/flirt -interp spline -in ${T1} -ref ${T1} -applyisoxfm 1 -out ${downsampled}/T1_1mm.nii.gz
        ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${T1} -r ${downsampled}/T1_1mm.nii.gz --premat=$FSLDIR/etc/flirtsch/ident.mat -o ${downsampled}/T1_1mm.nii.gz
        ${FSLDIR}/bin/applywarp --rel --interp=nn -i ${T1_brain} -r ${downsampled}/T1_1mm.nii.gz --premat=$FSLDIR/etc/flirtsch/ident.mat -o ${downsampled}/T1_brain_1mm.nii.gz 
        ${FSLDIR}/bin/fslmaths ${downsampled}/T1_1mm.nii.gz -div $Mean -mul 150 -abs ${downsampled}/T1_1mm.nii.gz
    else
        echo ${i} T1 downsampled to 1mm
    fi

#### recon-all ####
    #### autorecon1 (no skull stripping) ####
    recon-all -i ${T1_1mm} -subjid ${i} -sd ${FS} -motioncor -talairach -nuintensitycor -normalization 
    
    #### generate brain mask ####
    mri_convert ${T1_brain_1mm} ${FS}/mri/brainmask.mgz --conform
    mri_em_register -mask ${FS}/mri/brainmask.mgz ${i}/mri/nu.mgz $FREESURFER_HOME/average/RB_all_withskull_2016-05-10.vc700.gca ${FS}/mri/transforms/talairach_with_skull.lta
    mri_watershed -T1 -brain_atlas $FREESURFER_HOME/average/RB_all_withskull_2016-05-10.vc700.gca ${FS}/mri/transforms/talairach_with_skull.lta ${FS}/mri/T1.mgz ${FS}/mri/brainmask.auto.mgz 
    cp ${FS}/mri/brainmask.auto.mgz ${FS}/mri/brainmask.mgz 
    
    #### autorecon2 (turn off smooth2, inflate2, curvstats, and segstats)
    recon-all -subjid ${i} -sd ${FS}  -autorecon2 -nosmooth2 -noinflate2 -nocurvstats -nosegstats 

    #### call HCP script for high res WM surface placement & fine tuning T2 to T1 ####
    ${HCP_scripts}/FreeSurferHiresWhite.sh ${i} ${FS} ${T1} ${T2}

    #### intermediate recon-all steps ####
    recon-all -subjid ${i} -sd ${FS} -smooth2 -inflate2 -curvstats -sphere -surfreg -jacobian_white -avgcurv -cortparc 

    #### call HCP script for high res pial surface placement ####
    ${HCP_scripts}/FreeSurferHiresPial.sh ${i} ${FS} ${T1} ${T2}

    #### final recon-all steps ####
    recon-all -subjid ${i} -sd ${FS} -surfvolume -parcstats -cortparc2 -parcstats2 -cortparc3 -parcstats3 -cortribbon -segstats -aparc2aseg -wmparc -balabels -label-exvivo-ec

done




