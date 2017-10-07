FSLDIR=/usr/local/fsl
FREESURFER_HOME=/Applications/freesurfer
#export SUBJECTS_DIR=/home/yoobin_kwak/preproc_BCS_001_002
HCP_scripts=Pipelines-master/FreeSurfer/scripts

for i in $@
do
    FS=${i}/Preprocess_Structure/Freesurfer
    if [ ! -d ${FS} ]
    then
        mkdir ${FS}
    fi

    T1=${i}/Preprocess_Structure/PreFreesurfer/BiasField/T1wRestoredImage.nii.gz 
    T1_brain=${i}/Preprocess_Structure/PreFreesurfer/BiasField/T1wRestoredBrainImage.nii.gz 
    T2=${i}/Preprocess_Structure/PreFreesurfer/BiasField/T2wRestoredImage.nii.gz 

    export SUBJECTS_DIR=${FS}

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
        ${FSLDIR}/bin/flirt -interp spline -in ${T1} -ref ${T1} -applyisoxfm 1 -out ${T1_1mm}
        ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${T1} -r ${T1_1mm} --premat=$FSLDIR/etc/flirtsch/ident.mat -o ${T1_1mm}
        ${FSLDIR}/bin/applywarp --rel --interp=nn -i ${T1_brain} -r ${T1_1mm} --premat=$FSLDIR/etc/flirtsch/ident.mat -o ${T1_brain_1mm} 
        ${FSLDIR}/bin/fslmaths ${T1_1mm} -div $Mean -mul 150 -abs ${T1_1mm}
    else
        echo ${i} T1 downsampled to 1mm
    fi

#### recon-all ####
    fs_normalized=${FS}/${i}/mri/T1.mgz
    if [ ! -f ${fs_normalized} ]
    then
        #### autorecon1 (no skull stripping) ####
        recon-all -i ${T1_1mm} -subjid ${i} -sd ${FS} -motioncor -talairach -nuintensitycor -normalization 
    else
        echo ${i} FS 1st stage done
    fi

    fs_brainmask=${FS}/${i}/mri/brainmask.mgz
    if [ ! -f ${fs_brainmask} ]
    then
        #### generate brain mask (generates talairach_with_skull.log at the direcory level of the script ####
        mri_convert ${T1_brain_1mm} ${FS}/${i}/mri/brainmask.mgz --conform
        mri_em_register -mask ${FS}/${i}/mri/brainmask.mgz ${FS}/${i}/mri/nu.mgz $FREESURFER_HOME/average/RB_all_withskull_2008-03-26.gca ${FS}/${i}/mri/transforms/talairach_with_skull.lta
        mri_watershed -T1 -brain_atlas $FREESURFER_HOME/average/RB_all_withskull_2008-03-26.gca ${FS}/${i}/mri/transforms/talairach_with_skull.lta ${FS}/${i}/mri/T1.mgz ${FS}/${i}/mri/brainmask.auto.mgz 
        cp ${FS}/${i}/mri/brainmask.auto.mgz ${FS}/${i}/mri/brainmask.mgz 
    else 
        echo ${i} FS brainmask made
    fi

    fs_white_lh=${FS}/${i}/surf/lh.white 
    if [ ! -f ${fs_white_lh} ]
    then
        #### autorecon2 (turn off smooth2, inflate2, curvstats, and segstats)
        recon-all -subjid ${i} -sd ${FS} -autorecon2 -nosmooth2 -noinflate2 -nocurvstats -nosegstats 
    else
        echo ${i} FS 2nd stage done
    fi

    fs_T1w2T2=${FS}/${i}/mri/T2w_hires.nii.gz
    if [ ! -f ${fs_T1w2T2} ]
    then
        #### call HCP script for high res WM surface placement & fine tuning T2 to T1 (output log at the directory level of the script) ####
        bash ${HCP_scripts}/edit_FreeSurferHiresWhite.sh ${i} ${FS} ${T1} ${T2}
    else
        echo ${i} high resolution white surface placement and bbregister done
    fi

    fs_cortparc_lh=${FS}/${i}/label/lh.aparc.annot
    if [ ! -f ${fs_cortparc_lh} ]
    then
        #### intermediate recon-all steps ####
        recon-all -subjid ${i} -sd ${FS} -smooth2 -inflate2 -curvstats -sphere -surfreg -jacobian_white -avgcurv -cortparc 
    else
        echo ${i} FS stage 3 done
    fi

#    fs_ribbon1_dir=${FS}/${i}/mri/ribbon.postT2.pass1
#    if [ ! -d ${fs_ribbon1_dir} ]
#    then
#        #### call HCP script for high res pial surface placement ####
#        bash ${HCP_scripts}/edit_FreeSurferHiresPial.sh ${i} ${FS} ${T1} ${T2}
#    else
#        echo ${i} high resolution pial surface placement done
#    fi
#
#    fs_parcstats_lh=${FS}/${i}/stats/lh.aparc.a2005s.stats
#    if [ ! -f ${fs_parcstats_lh} ]
#    then
#        #### final recon-all steps ####
#        recon-all -subjid ${i} -sd ${FS} -surfvolume -parcstats -cortparc2 -parcstats2 -cortparc3 -parcstats3 -cortribbon -segstats -aparc2aseg -wmparc -balabels       ## -label-exvivo-ec option not used in version 6 (https://github.com/nipy/nipype/issues/1782)
#    else
#        echo ${i} FS final stage done
#    fi
#
done


