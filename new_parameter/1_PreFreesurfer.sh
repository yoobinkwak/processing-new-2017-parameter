FSLDIR=/usr/share/fsl/5.0
HCPDIR=Pipelines-master
CARET7DIR=workbench/bin_linux64

for i in $@
do
    preproc_struct=${i}/Preprocess_Structure
    if [ ! -d ${preproc_struct} ]
    then
        mkdir ${preproc_struct}
    fi
    
    pre_FS=${preproc_struct}/PreFreesurfer
    if [ ! -d ${pre_FS} ]
    then
        mkdir ${pre_FS}
    fi

    T1w=${pre_FS}/T1w
    if [ ! -d ${T1w} ]
    then
        mkdir ${T1w}
    fi

    T2w=${pre_FS}/T2w
    if [ ! -d ${T2w} ]
    then
        mkdir ${T2w}
    fi

#### ACPC alignment ####
    T1_raw=${i}/T1*/20*
    T2_raw=${i}/T2*/20*
    Reference_ACPC=/${FSLDIR}/data/standard/MNI152_T1_1mm.nii.gz

    count_t1w=`find ${T1w} -type f | wc -l`
    if [[ ${count_t1w} -lt 8 ]]
    then
        ${FSLDIR}/bin/robustfov -i ${T1_raw} -m ${T1w}/roi2full.mat -r ${T1w}/robustroi.nii.gz
        ${FSLDIR}/bin/convert_xfm -omat ${T1w}/full2roi.mat -inverse ${T1w}/roi2full.mat
        ${FSLDIR}/bin/flirt -interp spline -in ${T1w}/robustroi.nii.gz -ref ${Reference_ACPC} -omat ${T1w}/roi2std.mat -out ${T1w}/acpc_final.nii.gz -searchrx -180 180 -searchry -180 180 -searchrz -180 180
        ${FSLDIR}/bin/convert_xfm -omat ${T1w}/full2std.mat -concat ${T1w}/roi2std.mat ${T1w}/full2roi.mat
        ${FSLDIR}/bin/aff2rigid ${T1w}/full2std.mat ${T1w}/acpc_alignmnet.mat
        ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${T1_raw} -r ${Reference_ACPC} --premat=${T1w}/acpc_alignmnet.mat -o ${T1w}/ACPC_aligned.nii.gz
    else
        echo ${i} T1 aligned
    fi

    count_t2w=`find ${T2w} -type f | wc -l`
    if [[ ${count_t2w} -lt 8 ]]
    then
        ${FSLDIR}/bin/robustfov -i ${T2_raw} -m ${T2w}/roi2full.mat -r ${T2w}/robustroi.nii.gz
        ${FSLDIR}/bin/convert_xfm -omat ${T2w}/full2roi.mat -inverse ${T2w}/roi2full.mat
        ${FSLDIR}/bin/flirt -interp spline -in ${T2w}/robustroi.nii.gz -ref ${Reference_ACPC} -omat ${T2w}/roi2std.mat -out ${T2w}/acpc_final.nii.gz -searchrx -180 180 -searchry -180 180 -searchrz -180 180
        ${FSLDIR}/bin/convert_xfm -omat ${T2w}/full2std.mat -concat ${T2w}/roi2std.mat ${T2w}/full2roi.mat
        ${FSLDIR}/bin/aff2rigid ${T2w}/full2std.mat ${T2w}/acpc_alignmnet.mat
        ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${T2_raw} -r ${Reference_ACPC} --premat=${T2w}/acpc_alignmnet.mat -o ${T2w}/ACPC_aligned.nii.gz
    else
        echo ${i} T2 aligned
    fi

#### Brain extraction ####  CHECK: used T1 MNI references for T2
    T1_ACPC=${T1w}/ACPC_aligned.nii.gz
    Reference_forBrain=${HCPDIR}/global/templates/MNI152_T1_0.8mm.nii.gz
    Reference_forBrainMask=${HCPDIR}/global/templates/MNI152_T1_0.8mm_brain_mask.nii.gz 
    Reference_forBrain_2mm=${HCPDIR}/global/templates/MNI152_T1_2mm.nii.gz 
    Reference_forBrain_2mmMask=${HCPDIR}/global/templates/MNI152_T1_2mm_brain_mask_dil.nii.gz
    FNIRTConfig=${FSLDIR}/etc/flirtsch/T1_2_MNI152_2mm.cnf
    T1_brain=${T1w}/T1_brain.nii.gz
    if [ ! -f ${T1_brain} ]
    then
        ${FSLDIR}/bin/flirt -interp spline -dof 12 -in ${T1_ACPC} -ref ${Reference_forBrain_2mm} -omat ${T1w}/forBrain_roughlin.mat -out ${T1w}/T1ACPC_to_MNI_roughlin.nii.gz -nosearch
        ${FSLDIR}/bin/fnirt --in=${T1_ACPC} --ref=${Reference_forBrain_2mm} --aff=${T1w}/forBrain_roughlin.mat --refmask=${Reference_forBrain_2mmMask} --fout=${T1w}/forBrain_str2standard.nii.gz --jout=${T1w}/forBrain_NonlinearRegJacobians.nii.gz --refout=${T1w}/forBrain_IntensityModulatedT1.nii.gz --iout=${T1w}/T1ACPC_to_MNI_nonlin.nii.gz --logout=${T1w}/forBrain_NonlinearReg.txt --intout=${T1w}/forBrain_NonlinearIntensities.nii.gz --cout=${T1w}/forBrain_NonlinearReg.nii.gz --config=${FNIRTConfig}
        ${FSLDIR}/bin/applywarp --rel --interp=spline --in=${T1_ACPC} --ref=${Reference_forBrain} -w ${T1w}/forBrain_str2standard.nii.gz --out=${T1w}/T1ACPC_to_MNI_nonlin.nii.gz
        ${FSLDIR}/bin/invwarp --ref=${Reference_forBrain_2mm} -w ${T1w}/forBrain_str2standard.nii.gz -o ${T1w}/forBrain_standard2str.nii.gz
        ${FSLDIR}/bin/applywarp --rel --interp=nn --in=${Reference_forBrainMask} --ref=${T1_ACPC} -w ${T1w}/forBrain_standard2str.nii.gz -o ${T1w}/forBrain_FinalMask.nii.gz
        ${FSLDIR}/bin/fslmaths ${T1_ACPC} -mas ${T1w}/forBrain_FinalMask.nii.gz ${T1w}/T1_brain.nii.gz
    else
        echo ${i} T1 brain extracted
    fi

    T2_ACPC=${T2w}/ACPC_aligned.nii.gz
    T2_brain=${T2w}/T2_brain.nii.gz
    if [ ! -f ${T2_brain} ]
    then
        ${FSLDIR}/bin/flirt -interp spline -dof 12 -in ${T2_ACPC} -ref ${Reference_forBrain_2mm} -omat ${T2w}/forBrain_roughlin.mat -out ${T2w}/T2ACPC_to_MNI_roughlin.nii.gz -nosearch
        ${FSLDIR}/bin/fnirt --in=${T2_ACPC} --ref=${Reference_forBrain_2mm} --aff=${T2w}/forBrain_roughlin.mat --refmask=${Reference_forBrain_2mmMask} --fout=${T2w}/forBrain_str2standard.nii.gz --jout=${T2w}/forBrain_NonlinearRegJacobians.nii.gz --refout=${T2w}/forBrain_IntensityModulatedT2.nii.gz --iout=${T2w}/T2ACPC_to_MNI_nonlin.nii.gz --logout=${T2w}/forBrain_NonlinearReg.txt --intout=${T2w}/forBrain_NonlinearIntensities.nii.gz --cout=${T2w}/forBrain_NonlinearReg.nii.gz 
        ${FSLDIR}/bin/applywarp --rel --interp=spline --in=${T2_ACPC} --ref=${Reference_forBrain} -w ${T2w}/forBrain_str2standard.nii.gz --out=${T2w}/T2ACPC_to_MNI_nonlin.nii.gz
        ${FSLDIR}/bin/invwarp --ref=${Reference_forBrain_2mm} -w ${T2w}/forBrain_str2standard.nii.gz -o ${T2w}/forBrain_standard2str.nii.gz
        ${FSLDIR}/bin/applywarp --rel --interp=nn --in=${Reference_forBrainMask} --ref=${T2_ACPC} -w ${T2w}/forBrain_standard2str.nii.gz -o ${T2w}/forBrain_FinalMask.nii.gz
        ${FSLDIR}/bin/fslmaths ${T2_ACPC} -mas ${T2w}/forBrain_FinalMask.nii.gz ${T2w}/T2_brain.nii.gz
    else
        echo ${i} T2 brain extracted
    fi


#### T2 to T1 registration ####
    cp -r ${T1_ACPC} ${T1w}/T1.nii.gz
    cp -r ${T2_ACPC} ${T2w}/T2.nii.gz
    T1=${T1w}/T1.nii.gz
    T2=${T2w}/T2.nii.gz

    T2toT1_dir=${pre_FS}/T2toT1_Register
    if [ ! -d ${T2toT1_dir} ]
    then
        mkdir ${T2toT1_dir}
    fi
    T2_to_T1=${T2toT1_dir}/T2w2T1w.nii.gz
    if [ ! -f ${T2_to_T1} ]
    then
        ${FSLDIR}/bin/epi_reg --epi=${T2_brain} --t1=${T1} --t1brain=${T1_brain} --out=${T2toT1_dir}/T2w2T1w
        ${FSLDIR}/bin/applywarp --rel --interp=spline --in=${T2} --ref=${T1} --premat=${T2toT1_dir}/T2w2T1w.mat --out=${T2toT1_dir}/T2w2T1w
        ${FSLDIR}/bin/fslmaths ${T2toT1_dir}/T2w2T1w.nii.gz -add 1 ${T2toT1_dir}/T2w2T1w.nii.gz
    else
        echo ${i} T2 registerd to T1
    fi
    #### create files as if redout distortion correction was performed ####
    if [ ! -f ${T1w}/T2_acpc_dc.nii.gz ] 
    then
        cp ${T1} ${T1w}/T1_acpc_dc.nii.gz
        cp ${T1_brain} ${T1w}/T1_acpc_dc_brain.nii.gz
        mkdir ${T1w}/xfms 
        ${FSLDIR}/bin/fslmerge -t ${T1w}/xfms/T1_dc ${T1} ${T1} ${T1}
        ${FSLDIR}/bin/fslmaths ${T1w}/xfms/T1_dc -mul 0 ${T1w}/xfms/T1_dc 
        cp ${T2toT1_dir}/T2w2T1w.nii.gz ${T1w}/T2_acpc_dc.nii.gz
        ${FSLDIR}/bin/convertwarp --relout --rel -r ${T1w}/T2_acpc_dc.nii.gz -w ${T1w}/xfms/T1_dc --postmat=${T2toT1_dir}/T2w2T1w.mat --out=${T1w}/xfms/T2_reg_dc 
    else
        echo ${i} files created as if readout distortion correction was performed
    fi

#### Bias field correction #### requires Connectiome Workbench "wb_command" (line 131)
    bias_dir=${pre_FS}/BiasField
    if [ ! -d ${bias_dir} ]
    then
        mkdir ${bias_dir}
    fi
    BiasFieldSmoothingSigma=5
    Factor=0.5
    output_bias=${bias_dir}/BiasField.nii.gz
    if [ ! -f ${output_bias} ]
    then
        ${FSLDIR}/bin/fslmaths ${T1} -mul ${T2} -abs -sqrt ${bias_dir}/T1wmulT2w.nii.gz -odt float
        ${FSLDIR}/bin/fslmaths ${bias_dir}/T1wmulT2w.nii.gz -mas ${T1_brain} ${bias_dir}/T1wmulT2w_brain.nii.gz
        meanbrainval=`${FSLDIR}/bin/fslstats ${bias_dir}/T1wmulT2w_brain.nii.gz -M`
        ${FSLDIR}/bin/fslmaths ${bias_dir}/T1wmulT2w_brain.nii.gz -div $meanbrainval ${bias_dir}/T1wmulT2w_brain_norm.nii.gz
        ${FSLDIR}/bin/fslmaths ${bias_dir}/T1wmulT2w_brain_norm.nii.gz -bin -s ${BiasFieldSmoothingSigma} ${bias_dir}/SmoothNorm_s${BiasFieldSmoothingSigma}.nii.gz
        ${FSLDIR}/bin/fslmaths ${bias_dir}/T1wmulT2w_brain_norm.nii.gz -s ${BiasFieldSmoothingSigma} -div ${bias_dir}/SmoothNorm_s${BiasFieldSmoothingSigma}.nii.gz ${bias_dir}/T1wmulT2w_brain_norm_s${BiasFieldSmoothingSigma}.nii.gz
        ${FSLDIR}/bin/fslmaths ${bias_dir}/T1wmulT2w_brain_norm.nii.gz -div ${bias_dir}/T1wmulT2w_brain_norm_s${BiasFieldSmoothingSigma}.nii.gz ${bias_dir}/T1wmulT2w_brain_norm_modulate.nii.gz
        STD=`${FSLDIR}/bin/fslstats ${bias_dir}/T1wmulT2w_brain_norm_modulate.nii.gz -S`
        echo $STD
        MEAN=`${FSLDIR}/bin/fslstats ${bias_dir}/T1wmulT2w_brain_norm_modulate.nii.gz -M`
        echo $MEAN
        Lower=`echo "$MEAN - ($STD * $Factor)" | bc -l`
        echo $Lower
        ${FSLDIR}/bin/fslmaths ${bias_dir}/T1wmulT2w_brain_norm_modulate.nii.gz -thr $Lower -bin -ero -mul 255 ${bias_dir}/T1wmulT2w_brain_norm_modulate_mask.nii.gz
        ${CARET7DIR}/wb_command -volume-remove-islands ${bias_dir}/T1wmulT2w_brain_norm_modulate_mask.nii.gz ${bias_dir}/T1wmulT2w_brain_norm_modulate_mask.nii.gz
        ${FSLDIR}/bin/fslmaths ${bias_dir}/T1wmulT2w_brain_norm.nii.gz -mas ${bias_dir}/T1wmulT2w_brain_norm_modulate_mask.nii.gz -dilall ${bias_dir}/bias_raw.nii.gz -odt float
        ${FSLDIR}/bin/fslmaths ${bias_dir}/bias_raw.nii.gz -s ${BiasFieldSmoothingSigma} ${bias_dir}/BiasField.nii.gz
        ${FSLDIR}/bin/fslmaths ${T1} -div ${bias_dir}/BiasField.nii.gz -mas ${T1_brain} ${bias_dir}/T1wRestoredBrainImage.nii.gz -odt float
        ${FSLDIR}/bin/fslmaths ${T1} -div ${bias_dir}/BiasField.nii.gz ${bias_dir}/T1wRestoredImage.nii.gz -odt float 
        ${FSLDIR}/bin/fslmaths ${T2} -div ${bias_dir}/BiasField.nii.gz -mas ${T1_brain} ${bias_dir}/T2wRestoredBrainImage.nii.gz -odt float
        ${FSLDIR}/bin/fslmaths ${T2} -div ${bias_dir}/BiasField.nii.gz ${bias_dir}/T2wRestoredImage.nii.gz -odt float 
    else
        echo ${i} bias field corrected
    fi

#### MNI registration ####
    final_vol_register=${pre_FS}/Final_Register
    if [ ! -d ${final_vol_register} ]
    then
        mkdir ${final_vol_register}
    fi
    T1_restored=${bias_dir}/T1wRestoredImage.nii.gz 
    T1_restored_brain=${bias_dir}/T1wRestoredBrainImage.nii.gz 
    T2_restored=${bias_dir}/T2wRestoredImage.nii.gz 
    T2_restored_brain=${bias_dir}/T2wRestoredBrainImage.nii.gz 
    Reference_forBrainBrain=${HCPDIR}/global/templates/MNI152_T1_0.8mm_brain.nii.gz
    T1_final_reg=${final_vol_register}/T1wRestored_to_MNI_fnirt.nii.gz
    if [ ! -f ${T1_final_reg} ]
    then
        ${FSLDIR}/bin/flirt -interp spline -dof 12 -in ${T1_restored_brain} -ref ${Reference_forBrainBrain} -omat ${final_vol_register}/acpc2MNILinear.mat -out ${final_vol_register}/T1wRestored_Brain_to_MNI_flirt
        ${FSLDIR}/bin/fnirt --in=${T1_restored} --ref=${Reference_forBrain_2mm} --aff=${final_vol_register}/acpc2MNILinear.mat --refmask=${Reference_forBrain_2mmMask} --fout=${final_vol_register}/NonlinearTransform --jout=${final_vol_register}/NonlinearRegJacobians.nii.gz --refout=${final_vol_register}/IntensityModulatedT1.nii.gz --iout=${final_vol_register}/2mmReg.nii.gz --logout=${final_vol_register}/NonlinearReg.txt --intout=${final_vol_register}/NonlinearIntensities.nii.gz --cout=${final_vol_register}/NonlinearReg.nii.gz --config=${FNIRTConfig}
        ${FSLDIR}/bin/invwarp -w ${final_vol_register}/NonlinearTransform -o ${final_vol_register}/InvNonlinearTransform -r ${Reference_forBrain_2mm}
        ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${T1} -r ${Reference_forBrain} -w ${final_vol_register}/NonlinearTransform -o ${final_vol_register}/T1w_to_MNI_fnirt
        ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${T1_restored} -r ${Reference_forBrain} -w ${final_vol_register}/NonlinearTransform -o ${final_vol_register}/T1wRestored_to_MNI_fnirt
        ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${T1_restored_brain} -r ${Reference_forBrain} -w ${final_vol_register}/NonlinearTransform -o ${final_vol_register}/T1wRestoredBrain_to_MNI_fnirt
        ${FSLDIR}/bin/fslmaths ${final_vol_register}/T1wRestored_to_MNI_fnirt -mas ${final_vol_register}/T1wRestoredBrain_to_MNI_fnirt ${final_vol_register}/T1wRestoredBrain_to_MNI_fnirt
        ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${T2} -r ${Reference_forBrain} -w ${final_vol_register}/NonlinearTransform -o ${final_vol_register}/T2w_to_MNI_fnirt
        ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${T2_restored} -r ${Reference_forBrain} -w ${final_vol_register}/NonlinearTransform -o ${final_vol_register}/T2wRestored_to_MNI_fnirt
        ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${T2_restored_brain} -r ${Reference_forBrain} -w ${final_vol_register}/NonlinearTransform -o ${final_vol_register}/T2wRestoredBrain_to_MNI_fnirt
        ${FSLDIR}/bin/fslmaths ${final_vol_register}/T2wRestored_to_MNI_fnirt -mas ${final_vol_register}/T2wRestoredBrain_to_MNI_fnirt ${final_vol_register}/T2wRestoredBrain_to_MNI_fnirt
    else
        echo ${i} T1 final registered
    fi


done
