## Set directory paths
FSLDIR=/usr/local/fsl
HCPPIPEDIR=/Users/yoobin_kwak/preproc_BCS_001_002/Pipelines-master
HCPPIPEDIR_FS=${HCPPIPEDIR}/FreeSurfer/scripts
PipelineScripts=${HCPPIPEDIR_FS}
FREESURFER_HOME=/Applications/freesurfer

for i in $@
do
    SubjectID=${i}
    SubjectDIR=${i}/Preprocess_Structure/Freesurfer
    if [ ! -d ${SubjectDIR} ]
    then
        mkdir ${SubjectDIR}
    fi
    T1wImage=${i}/Preprocess_Structure/PreFreesurfer/T1w/T1w_acpc_dc_restore.nii.gz
    T1wImageBrain=${i}/Preprocess_Structure/PreFreesurfer/T1w/T1w_acpc_dc_restore_brain.nii.gz
    T2wImage=${i}/Preprocess_Structure/PreFreesurfer/T1w/T2w_acpc_dc_restore.nii.gz

    #여기 아래 코드 바꿀 수 있음
    T1wImageFile=`remove_ext $T1wImage`;
    T1wImageBrainFile=`remove_ext $T1wImageBrain`;

    export SUBJECTS_DIR=${SubjectDIR}

#### Make Spline Interpolated Downsample to 1mm ####
    if [ ! -f "$T1wImageFile"_1mm.nii.gz ] ; then
        echo "${i} FREESURFER : Making Spline Interpolated Downsample to 1mm"
        Mean=`fslstats $T1wImageBrain -M`
        flirt -interp spline -in "$T1wImage" -ref "$T1wImage" -applyisoxfm 1 -out "$T1wImageFile"_1mm.nii.gz
        applywarp --rel --interp=spline -i "$T1wImage" -r "$T1wImageFile"_1mm.nii.gz --premat=$FSLDIR/etc/flirtsch/ident.mat -o "$T1wImageFile"_1mm.nii.gz
        applywarp --rel --interp=nn -i "$T1wImageBrain" -r "$T1wImageFile"_1mm.nii.gz --premat=$FSLDIR/etc/flirtsch/ident.mat -o "$T1wImageBrainFile"_1mm.nii.gz
        fslmaths "$T1wImageFile"_1mm.nii.gz -div $Mean -mul 150 -abs "$T1wImageFile"_1mm.nii.gz
    else
        echo "${i} FREESURFER : Input T1w Downsampled to 1mm"
    fi

#### Initial Recon-all Steps ####
    if [ ! -f ${SubjectDIR}/${i}/mri/T1.mgz ] ; then         ## normalization output
        echo "${i} FREESURFER : Runnung autorecon1 steps with the exception of -skullstrip" 
        # Call recon-all with flags that are part of "-autorecon1", with the exception of -skullstrip. 
        # -skullstrip of FreeSurfer not reliable for Phase II data because of poor FreeSurfer mri_em_register registrations with Skull on, so run registration with PreFreeSurfer masked data and then generate brain mask as usual.
        recon-all -i "$T1wImageFile"_1mm.nii.gz -subjid $SubjectID -sd $SubjectDIR -motioncor -talairach -nuintensitycor -normalization 
    else
        echo "${i} FREESURFER : 1st (out of 4) recon-all Stage Completed  "
    fi

#### Generate brain mask ####
    if [ ! -f ${SubjectDIR}/${i}/mri/brainmask.mgz ] ; then
        echo "${i} FREESURFER : Generating Brain Mask "
        mri_convert "$T1wImageBrainFile"_1mm.nii.gz "$SubjectDIR"/"$SubjectID"/mri/brainmask.mgz --conform
        mri_em_register -mask "$SubjectDIR"/"$SubjectID"/mri/brainmask.mgz "$SubjectDIR"/"$SubjectID"/mri/nu.mgz $FREESURFER_HOME/average/RB_all_2008-03-26.gca "$SubjectDIR"/"$SubjectID"/mri/transforms/talairach_with_skull.lta
        mri_watershed -T1 -brain_atlas $FREESURFER_HOME/average/RB_all_withskull_2008-03-26.gca "$SubjectDIR"/"$SubjectID"/mri/transforms/talairach_with_skull.lta "$SubjectDIR"/"$SubjectID"/mri/T1.mgz "$SubjectDIR"/"$SubjectID"/mri/brainmask.auto.mgz
        cp "$SubjectDIR"/"$SubjectID"/mri/brainmask.auto.mgz "$SubjectDIR"/"$SubjectID"/mri/brainmask.mgz
        mv talairach_with_skull.log ${SubjectDIR}/${i}/
    else 
        echo "${i} FREESURFER : Brain Mask Generated"
    fi

#### Call recon-all to run most of the "-autorecon2" stages, but turning off smooth2, inflate2, curvstats, and segstats stages ####
    if [ ! -f ${SubjectDIR}/${i}/surf/lh.white ] ; then         ## finalsurfs output (check if this is the last output generated from this stage)
        echo "${i} FREESURFER : Runnung autorecon2 steps with few exceptions"
        recon-all -subjid $SubjectID -sd $SubjectDIR -autorecon2 -nosmooth2 -noinflate2 -nocurvstats -nosegstats
    else
        echo "${i} FREESURFER : 2nd (out of 4) recon-all Stage Completed  "
    fi

#### Highres white stuff and Fine Tune T2w to T1w Reg ####
    if [ ! -f ${SubjectDIR}/${i}/mri/T2w_hires.nii.gz ] ; then
        echo "${i} : FREESURFER Processing High Resolution WM and Fine Tuning T2w to T1w Registration"
        bash ${PipelineScripts}/edit_FreeSurferHiresWhite.sh "$SubjectID" "$SubjectDIR" "$T1wImage" "$T2wImage"
        mv *h.white.deformed.out ${SubjectDIR}/${i}/
    else
        echo "${i} FREESURFER : High Resolution White Surface Placement and bbregister Completed"
    fi

#### Intermediate Recon-all Steps ####
    if [ ! -f ${SubjectDIR}/${i}/label/lh.aparc.annot ] ; then      ## cortpar output
        echo "${i} FREESURFER : Runnung intermediate recon-all steps"
        recon-all -subjid $SubjectID -sd $SubjectDIR -smooth2 -inflate2 -curvstats -sphere -surfreg -jacobian_white -avgcurv -cortparc
    else
        echo "${i} FREESURFER : 3rd (out of 4) recon-all Stage  Completed"
    fi

#### Highres pial stuff (this module adjusts the pial surface based on the the T2w image) ####
    if [ ! -f ${SubjectDIR}/${i}/surf/rh.thickness.postT2.pass2 ] ; then        
        echo "${i} FREESURFER : Running High Resolution Pial Surface"

        bash ${PipelineScripts}/edit_FreeSurferHiresPial.sh "$SubjectID" "$SubjectDIR" "$T1wImage" "$T2wImage"
    else
        echo "${i} FREESURFER : High Resolution Pial Surface Placement Completed"
    fi

#### Final Recon-all Steps ####
    if [ ! -f ${SubjectDIR}/${i}/stats/lh.entorhinal_exvivo.stats ] ; then
        echo "${i} FREESURFER : Running final recon-all steps"

        recon-all -subjid $SubjectID -sd $SubjectDIR -surfvolume -parcstats -cortparc2 -parcstats2 -cortparc3 -parcstats3 -cortribbon -segstats -aparc2aseg -wmparc -balabels -label-exvivo-ec
    else
        echo "${i} FREESURFER : 4th (out of 4) recon-all Stage  Completed"
    fi

done


