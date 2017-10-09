FSLDIR=/usr/share/fsl/5.0 
HCPDIR=Pipelines-master
HCPPIPEDIR_Global=${HCPDIR}/global/scripts

for i in $@
do
    preproc_func=${i}/Preprocess_Function
    if [ ! -d ${preproc_func} ]
    then
        mkdir ${preproc_func}
    fi

    REST_MB4_LR=${i}/REST_MB4_LR_SBREF_0*/20* 
    REST_SBREF_LR=${i}/REST_MB4_LR_SBREF_SBREF*/20*

#### Motion correction (trying 2 methods) ####    
    mc_dir1=${preproc_func}/MotionCorrection_mcflirt
    if [ ! -d ${mc_dir1} ]
    then
        mkdir ${mc_dir1}
    fi
    mc_matrix_dir1=${mc_dir1}/MotionMatrices
    if [ ! -d ${mc_matrix_dir1} ]
    then
        mkdir ${mc_matrix_dir1}
    fi
    if [ ! -f ${mc_dir1}/motion_corrected.nii.gz ]
    then
        bash ${HCPDIR}/fMRIVolume/scripts/edit_MotionCorrection.sh ${mc_dir1} ${REST_MB4_LR} ${REST_SBREF_LR} motion_corrected motion_regressors ${mc_matrix_dir1} motion_matrix MCFLIRT
    else
        echo ${i} motion corrected with mcflirt
    fi

    mc_dir2=${preproc_func}/MotionCorrection_flirt
    if [ ! -d ${mc_dir2} ]
    then
        mkdir ${mc_dir2}
    fi
    mc_matrix_dir2=${mc_dir2}/MotionMatrices
    if [ ! -d ${mc_matrix_dir2} ]
    then
        mkdir ${mc_matrix_dir2}
    fi
    if [ ! -f ${mc_dir2}/motion_corrected.nii.gz ]
    then
        bash ${HCPDIR}/fMRIVolume/scripts/edit_MotionCorrection.sh ${mc_dir2} ${REST_MB4_LR} ${REST_SBREF_LR} motion_corrected motion_regressors ${mc_matrix_dir2} motion_matrix FLIRT
    else
        echo ${i} motion corrected with flirt
    fi

#### EPI Distortion Correction (from TopupPrerpocessingAll.sh) ####
    distortion=${preproc_func}/EPIDistortionCorrection
    if [ ! -d ${distortion} ]
    then
        mkdir ${distortion}
    fi
    SE_Pos=${i}/REST_MB1_BLIP_RL*/20*       ## "Pos" bc in PE dir (i.e., RL); RL/X/1
    SE_Neg=${i}/REST_MB1_BLIP_LR*/20*        ## "Neg" bc in opposite PE dir (i.e., LR); LR/X-/-1
    echo_spacing=0.00069        ## in sec
    UnWarpDir=x
    TopupConfig=${HCPDIR}/global/config/b02b0.cnf
    source $HCPPIPEDIR_Global/log.shlib
    
    #### check dimensions of phase versus sbref images ####
    if [[ `fslhd ${SE_Neg} | grep '^dim[123]'` != `fslhd ${REST_SBREF_LR} | grep '^dim[123]'` ]]
    then
        log_Msg "Error: Spin echo fieldmap has different dimensions than scout image, this requires a manual fix"
    fi
    #### check that the spin echo images match ####
    if [[ `fslhd ${SE_Neg} | grep '^dim[123]'` != `fslhd ${SE_Pos} | grep '^dim[123]'` ]]
    then
        log_Msg "Error: Spin echo fieldmap images have different dimensions!"
    fi
    #### raw input files and process file names to "match"  gradient-distortion-corrected files (i.e., these files are not corrected for gradient distrotion) ####
    if [ ! -f ${distortion}/PhaseOne.nii.gz ]
    then
        ${FSLDIR}/bin/imcp ${SE_Neg} ${distortion}/PhaseOne.nii.gz
        ${FSLDIR}/bin/imcp ${SE_Pos} ${distortion}/PhaseTwo.nii.gz
        ${FSLDIR}/bin/imcp ${REST_SBREF_LR} ${distortion}/SBRef.nii.gz
    else
        echo ${i} topup input files exist
    fi
    if [ ! -f ${distortion}/PhaseOne_gdc.nii.gz ] 
    then
        cp ${distortion}/PhaseOne.nii.gz ${distortion}/PhaseOne_gdc.nii.gz
        cp ${distortion}/PhaseTwo.nii.gz ${distortion}/PhaseTwo_gdc.nii.gz
        fslmerge -t ${distortion}/BothPhases ${distortion}/PhaseOne_gdc ${distortion}/PhaseTwo_gdc
        fslmaths ${distortion}/PhaseOne_gdc.nii.gz -mul 0 -add 1 ${distortion}/Mask
    else
        echo ${i} file names processed to match gradient-distortion-corrected file names
    fi
    
    #### Run Topup to estimate distortion in the phase encoding direction ####
    dimtOne=`${FSLDIR}/bin/fslval ${distortion}/PhaseOne dim4`
    dimtTwo=`${FSLDIR}/bin/fslval ${distortion}/PhaseTwo dim4`
    if [ ! -f ${distortion}/Magnitudes.nii.gz ]
    then
        #### Set up acqparams.txt for Topup ####
        txtfname=${distortion}/acqparams.txt
        if [ ! -f $txtfname ]
        then
            dimx=`${FSLDIR}/bin/fslval ${distortion}/PhaseOne dim1`
            nPEsteps=$(($dimx - 1))
            ro_time=`echo "scale=6; ${echo_spacing} * ${nPEsteps}" | bc -l`        ## check if grappa factor has to be included (current version does not include grappa factor)
            log_Msg "Total readout time is $ro_time secs"
            j=1
            while [ $j -le $dimtOne ]
            do
                echo "-1 0 0 $ro_time" >> $txtfname
                ShiftOne="x-"
                j=`echo "$j + 1" | bc`
            done
            j=1
            while [ $j -le $dimtTwo ] 
            do
                echo "1 0 0 $ro_time" >> $txtfname
                ShiftTwo="x"
                j=`echo "$j + 1" | bc`
            done
        else
            echo ${i} acqparams.txt created fro Topup
        fi
        #### Extrapolate the existing values beyond the mask (adding 1 just to avoid smoothing inside the mask) ####
        ${FSLDIR}/bin/fslmaths ${distortion}/BothPhases -abs -add 1 -mas ${distortion}/Mask -dilM -dilM -dilM -dilM -dilM ${distortion}/BothPhases
        #### Run Topup ####
        ${FSLDIR}/bin/topup --imain=${distortion}/BothPhases --datain=$txtfname --config=${TopupConfig} --out=${distortion}/Coefficents --iout=${distortion}/Magnitudes --fout=${distortion}/TopupField --dfout=${distortion}/WarpField --rbmout=${distortion}/MotionMatrix --jacout=${distortion}/Jacobian -v
    else
        echo ${i} Topup for REST completed
    fi

    #### Undistort #### 
    if [ ! -f ${distortion}/Magnitude_brain.nii.gz ]
    then
        if [[ $UnWarpDir = "x" ]]
        then
            #### select the first volume from PhaseTwo #### 
            VolumeNumber=$(($dimtOne + 1))
            vnum=`${FSLDIR}/bin/zeropad $VolumeNumber 2`
            #### register scout to SE input (PhaseTwo) + combine motion and distortion correction ####
            ${FSLDIR}/bin/flirt -dof 6 -interp spline -in ${distortion}/SBRef.nii.gz -ref ${distortion}/PhaseTwo_gdc -omat ${distortion}/SBRef2PhaseTwo_gdc.mat -out ${distortion}/SBRef2PhaseTwo_gdc
            ${FSLDIR}/bin/convert_xfm -omat ${distortion}/SBRef2WarpField.mat -concat ${distortion}/MotionMatrix_${vnum}.mat ${distortion}/SBRef2PhaseTwo_gdc.mat
            ${FSLDIR}/bin/convertwarp --relout --rel -r ${distortion}/PhaseTwo_gdc --premat=${distortion}/SBRef2WarpField.mat --warp1=${distortion}/WarpField_${vnum} --out=${distortion}/WarpField.nii.gz
            ${FSLDIR}/bin/imcp ${distortion}/Jacobian_${vnum}.nii.gz ${distortion}/Jacobian.nii.gz
            SBRefPhase=Two
        elif [[ $UnWarpDir = "x-" || $UnWarpDir = "-x" ]]
        then
            #### select the first volume from PhaseOne #### 
            VolumeNumber=$((0 + 1))
            vnum=`${FSLDIR}/bin/zeropad $VolumeNumber 2`
            #### register scout to SE input (PhaseOne) + combine motion and distortion correction ####
            ${FSLDIR}/bin/flirt -dof 6 -interp spline -in ${distortion}/SBRef.nii.gz -ref ${distortion}/PhaseOne_gdc -omat ${distortion}/SBRef2PhaseOne_gdc.mat -out ${distortion}/SBRef2PhaseOne_gdc
            ${FSLDIR}/bin/convert_xfm -omat ${distortion}/SBRef2WarpField.mat -concat ${distortion}/MotionMatrix_${vnum}.mat ${distortion}/SBRef2PhaseOne_gdc.mat
            ${FSLDIR}/bin/convertwarp --relout --rel -r ${distortion}/PhaseOne_gdc --premat=${distortion}/SBRef2WarpField.mat --warp1=${distortion}/WarpField_${vnum} --out=${distortion}/WarpField.nii.gz
            ${FSLDIR}/bin/imcp ${distortion}/Jacobian_${vnum}.nii.gz ${distortion}/Jacobian.nii.gz
            SBRefPhase=One
        fi
        #### PhaseTwo (first vol) - warp and Jacobian modulate to get distortion corrected output ####
        VolumeNumber=$(($dimtOne + 1))
        vnum=`${FSLDIR}/bin/zeropad $VolumeNumber 2`
        ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${distortion}/PhaseTwo_gdc -r ${distortion}/PhaseTwo_gdc --premat=${distortion}/MotionMatrix_${vnum}.mat -w ${distortion}/WarpField_${vnum} -o ${distortion}/PhaseTwo_gdc_dc
        ${FSLDIR}/bin/fslmaths ${distortion}/PhaseTwo_gdc_dc -mul ${distortion}/Jacobian_${vnum} ${distortion}/PhaseTwo_gdc_dc_jac
        #### PhaseOne (first vol) - warp and Jacobian modulate to get distortion corrected output ####
        VolumeNumber=$((0 + 1))
        vnum=`${FSLDIR}/bin/zeropad $VolumeNumber 2`
        ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${distortion}/PhaseOne_gdc -r ${distortion}/PhaseOne_gdc --premat=${distortion}/MotionMatrix_${vnum}.mat -w ${distortion}/WarpField_${vnum} -o ${distortion}/PhaseOne_gdc_dc
        ${FSLDIR}/bin/fslmaths ${distortion}/PhaseOne_gdc_dc -mul ${distortion}/Jacobian_${vnum} ${distortion}/PhaseOne_gdc_dc_jac
        #### Scout - warp and Jacobian modulate to get distortion corrected output ####
        ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${distortion}/SBRef.nii.gz -r ${distortion}/SBRef.nii.gz -w ${distortion}/WarpField.nii.gz -o ${distortion}/SBRef_dc.nii.gz
        ${FSLDIR}/bin/fslmaths ${distortion}/SBRef_dc.nii.gz -mul ${distortion}/Jacobian.nii.gz ${distortion}/SBRef_dc_jac.nii.gz
        #### Calculate Equivalent Field Map ####
#        ${FSLDIR}/bin/fslmaths ${distortion}/TopupField -mul 6.283 ${distortion}/TopupField
        ${FSLDIR}/bin/fslmaths ${distortion}/Magnitudes.nii.gz -Tmean ${distortion}/Magnitude.nii.gz
        ${FSLDIR}/bin/bet ${distortion}/Magnitude ${distortion}/Magnitude_brain -f 0.35 -m #Brain extract the magnitude image
    else
        echo ${i} REST corrected for distortion on the phase encoding direction
    fi












    














done


