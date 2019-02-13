FSLDIR=/usr/local/fsl
HCPPIPEDIR=/Users/yoobin_kwak/preproc_BCS_001_002/Pipelines-master
HCPPIPEDIR_fMRIVol=${HCPPIPEDIR}/fMRIVolume/scripts
PipelineScripts=${HCPPIPEDIR_fMRIVol}
HCPPIPEDIR_Global=${HCPPIPEDIR}/global/scripts
GlobalScripts=${HCPPIPEDIR_Global}
FREESURFER_HOME=/Applications/freesurfer
CARET7DIR=/Users/yoobin_kwak/preproc_BCS_001_002/workbench/bin_macosx64

for i in $@ ; do
    DwellTime=0.00069        ## in sec
    grappa_factor=2
    UnwarpDir=x-         
    FinalfMRIResolution=2 
    DistortionCorrection=TOPUP         
    BiasCorrection=LEGACY        #SEBASED (use bias field derived from spin echo, must also use TOPUP, and require PostFreesurfer), LEGACY (use the bias field derived from T1w and T2w images) or NONE         
    GradientDistortionCoeffs=NONE
    TopupConfig=${HCPPIPEDIR}/global/config/b02b0.cnf
    dof=6       ## Default
    UseJacobian=true            #true or false (NOTE: the jacobian option only applies the jacobian of the distortion corrections to the fMRI data, and NOT from the nonlinear T1 to template registration)
    #MotionCorrectionType=      #MCFLIRT (default) or FLIRT

    #Subject=${i}
    Subject=${i}/baseline
    outLoc=${Subject}/preprocessed
    #fMRITimeSeries=${i}/REST_MB4_LR_SBREF_0*/20* 
    fMRITimeSeries=${outLoc}/REST/*gz
    #fMRIScout=${i}/REST_MB4_LR_SBREF_SBREF*/20*
    fMRIScout=${outLoc}/REST_LR_SBREF/*gz
    #SpinEchoPhaseEncodeNegative=${i}/REST_MB1_BLIP_LR*/20*        ## "Neg" bc in opposite PE dir (i.e., LR); LR/X-/-1
    #SpinEchoPhaseEncodePositive=${i}/REST_MB1_BLIP_RL*/20*       ## "Pos" bc in PE dir (i.e., RL); RL/X/1
    SpinEchoPhaseEncodeNegative=${outLoc}/REST_LR/*gz       ## "Neg" bc in opposite PE dir (i.e., LR); LR/X-/-1
    SpinEchoPhaseEncodePositive=${outLoc}/REST_RL/*gz       ## "Pos" bc in PE dir (i.e., RL); RL/X/1
    #NameOffMRI=Preprocess_Function
    NameOffMRI=Funct_preproc
    NamingOffMRI=REST       ## I added (take caution in the proceeding naming conventions)
    T1wImage="T1w_acpc_dc"
    T1wRestoreImage="T1w_acpc_dc_restore"
    T1wRestoreImageBrain="T1w_acpc_dc_restore_brain"
    ResultsFolder="Results"
    BiasField="BiasField_acpc_dc"
    BiasFieldMNI="BiasField"
    T1wAtlasName="T1w_restore"
    MovementRegressor="Movement_Regressors" #No extension, .txt appended
    MotionMatrixFolder="MotionMatrices"
    MotionMatrixPrefix="MAT_"
    FieldMapOutputName="FieldMap"
    MagnitudeOutputName="Magnitude"
    MagnitudeBrainOutputName="Magnitude_brain"
    ScoutName="Scout"
    OrigScoutName="${ScoutName}_orig"
    OrigTCSName="${NamingOffMRI}_orig"
    FreeSurferBrainMask="brainmask_fs"
    fMRI2strOutputTransform="${NamingOffMRI}2str"
    RegOutput="Scout2T1w"
    AtlasTransform="acpc_dc2standard"
    OutputfMRI2StandardTransform="${NamingOffMRI}2standard"
    Standard2OutputfMRITransform="standard2${NamingOffMRI}"
    QAImage="T1wMulEPI"
    JacobianOut="Jacobian"
    SubjectFolder=$Subject
    sebasedBiasFieldMNI="$AtlasSpaceFolder/Results/$NamingOffMRI/${NamingOffMRI}_sebased_bias.nii.gz"    #note, this file doesn't exist yet, gets created by ComputeSpinEchoBiasField.sh during DistortionCorrectionAnd...
    #T1wFolder=${i}/Preprocess_Structure/PreFreesurfer/T1w
    T1wFolder=${outLoc}/Struct_preproc/T1w
    #AtlasSpaceFolder=${i}/Preprocess_Structure/PreFreesurfer/MNINonLinear
    AtlasSpaceFolder=${outLoc}/Struct_preproc/MNINonLinear
    if [ ! -e $AtlasSpaceFolder/$ResultsFolder/$NamingOffMRI ] ; then
        mkdir -p $AtlasSpaceFolder/$ResultsFolder/$NamingOffMRI
    fi
    if [ ! -e ${T1wFolder}/Results/$NamingOffMRI ] ; then
        mkdir -p ${T1wFolder}/Results/$NamingOffMRI
    fi
    #fMRIFolder=$Subject/$NameOffMRI
    fMRIFolder=${outLoc}/$NameOffMRI
    if [ ! -e $fMRIFolder ] ; then
        mkdir $fMRIFolder
    fi
    if [ ! -f $fMRIFolder/"$OrigTCSName".nii.gz ] ; then
        cp $fMRITimeSeries $fMRIFolder/"$OrigTCSName".nii.gz
    fi
    if [ ! -f $fMRIFolder/"$OrigScoutName".nii.gz ] ; then
        cp $fMRIScout $fMRIFolder/"$OrigScoutName".nii.gz
    fi

#### Gradient Distortion Correction of fMRI (from GenericfMRIVolumeProcessingPipeline.sh) ####
    if [ $GradientDistortionCoeffs = "NONE" ] ; then
        echo "${i} Resting Functional Preprocessing : NOT PERFORMING GRADIENT DISTORTION CORRECTION"
        if [ ! -f $fMRIFolder/"$NamingOffMRI"_gdc_warp_jacobian.nii.gz ] ; then
            ${FSLDIR}/bin/imcp $fMRIFolder/$OrigTCSName $fMRIFolder/"$NamingOffMRI"_gdc
            ${FSLDIR}/bin/fslroi $fMRIFolder/"$NamingOffMRI"_gdc $fMRIFolder/"$NamingOffMRI"_gdc_warp 0 3
            ${FSLDIR}/bin/fslmaths $fMRIFolder/"$NamingOffMRI"_gdc_warp -mul 0 $fMRIFolder/"$NamingOffMRI"_gdc_warp
            ${FSLDIR}/bin/imcp $fMRIFolder/$OrigScoutName $fMRIFolder/"$ScoutName"_gdc
            #make fake jacobians of all 1s, for completeness
            ${FSLDIR}/bin/fslmaths $fMRIFolder/$OrigScoutName -mul 0 -add 1 $fMRIFolder/"$ScoutName"_gdc_warp_jacobian
            ${FSLDIR}/bin/fslroi $fMRIFolder/"$NamingOffMRI"_gdc_warp $fMRIFolder/"$NamingOffMRI"_gdc_warp_jacobian 0 1
            ${FSLDIR}/bin/fslmaths $fMRIFolder/"$NamingOffMRI"_gdc_warp_jacobian -mul 0 -add 1 $fMRIFolder/"$NamingOffMRI"_gdc_warp_jacobian
        elif [  -e $fMRIFolder/"$NamingOffMRI"_gdc_warp_jacobian.nii.gz ] ; then
            correct_outputs=5
            check_outputs=`ls $fMRIFolder/*_gdc*nii.gz | wc -l`
            if [ $check_outputs -eq $correct_outputs ] ; then
                echo "      Scout and TimeSeries Images Processed as if Gradient Distortion Correction and Jacobian were Applied "
            else
                echo "      ERROR Porcessing Scout and TimeSeries Images as if Gradient Distortion Correction and Jacobian were Applied "
            fi
        fi
    else
        echo "${i} Resting Functional Preprocess : Performing Gradient Distortion Correction on Scout and TimeSeries Images"
        if [ ! -d $fMRIFolder/GradientDistortionUnwarp ] || [ ! -d $fMRIFolder/"$ScoutName"_GradientDistortionUnwarp ] ; then
            mkdir $fMRIFolder/GradientDistortionUnwarp
            bash $GlobalScripts/GradientDistortionUnwarp.sh --workingdir=$fMRIFolder/GradientDistortionUnwarp --coeffs=$GradientDistortionCoeffs --in=$fMRIFolder/$OrigTCSName --out=$fMRIFolder/"$NamingOffMRI"_gdc --owarp=$fMRIFolder/"$NamingOffMRI"_gdc_warp
            mkdir  $fMRIFolder/"$ScoutName"_GradientDistortionUnwarp
            bash $GlobalScripts/GradientDistortionUnwarp.sh --workingdir=$fMRIFolder/"$ScoutName"_GradientDistortionUnwarp --coeffs=$GradientDistortionCoeffs --in=$fMRIFolder/$OrigScoutName --out=$fMRIFolder/"$ScoutName"_gdc --owarp=$fMRIFolder/"$ScoutName"_gdc_warp
            if [[ $UseJacobian == "true" ]] ; then
                ${FSLDIR}/bin/fslmaths $fMRIFolder/"$NamingOffMRI"_gdc -mul $fMRIFolder/"$NamingOffMRI"_gdc_warp_jacobian $fMRIFolder/"$NamingOffMRI"_gdc
                ${FSLDIR}/bin/fslmaths $fMRIFolder/"$ScoutName"_gdc -mul $fMRIFolder/"$ScoutName"_gdc_warp_jacobian $fMRIFolder/"$ScoutName"_gdc
            fi
        fi
    fi

#### Motion Correction ####
    MC_dir1=$fMRIFolder/MotionCorrection_MCFLIRT
    if [ ! -e ${MC_dir1} ]; then
        mkdir ${MC_dir1}
    fi
#    MC_dir2=$fMRIFolder/MotionCorrection_FLIRT
#    if [ ! -e ${MC_dir2} ]; then
#        mkdir ${MC_dir2}
#    fi
    InputfMRI=$fMRIFolder/"$NamingOffMRI"_gdc 
    Scout=$fMRIFolder/"$ScoutName"_gdc

    if [ ! -f ${MC_dir1}/"$NamingOffMRI"_mc_mcflirt.nii.gz ] ; then
        echo "${i} Resting Functional Preprocess : Performing motion correction with MCFLIRT"
        bash ${HCPPIPEDIR_fMRIVol}/edit_MotionCorrection.sh ${MC_dir1} ${InputfMRI} ${Scout} "$NamingOffMRI"_mc_mcflirt "$MovementRegressor"_mcflirt $fMRIFolder/"$MotionMatrixFolder"_mcflirt $MotionMatrixPrefix MCFLIRT
    else
        echo "${i} Resting Functional Preprocess : Motion correction with MCFLIRT Completed"
    fi

#    if [ ! -f ${MC_dir2}/"$NamingOffMRI"_mc_flirt.nii.gz ] ; then
#        echo "${i} Resting Functional Preprocess : Performing motion correction with FLIRT"
#        bash ${HCPPIPEDIR_fMRIVol}/edit_MotionCorrection.sh ${MC_dir2} ${InputfMRI} ${Scout} "$NamingOffMRI"_mc_flirt "$MovementRegressor"_flirt $fMRIFolder/"$MotionMatrixFolder"_flirt $MotionMatrixPrefix FLIRT
#    else
#        echo "${i} Resting Functional Preprocess : Motion correction with FLIRT Completed"
#    fi

#### EPI Distortion Correction and EPI to T1w Registration (from DistortionCorrectionAndEPIToT1wReg_FLIRTBBRAndFreeSurferBBRbased.sh, lines 136~700) ####
    DCFolder=${fMRIFolder}/DistortionCorrectionAndEPIToT1wReg_FLIRTBBRAndFreeSurferBBRbased
    if [ ! -e ${DCFolder} ] ; then
        mkdir ${DCFolder}
    fi
    ScoutInputName=${fMRIFolder}/${ScoutName}_gdc 
    ScoutInputFile=`basename $ScoutInputName`
    T1wImage=${T1wFolder}/${T1wImage} 
    T1wRestoreImage=${T1wFolder}/${T1wRestoreImage} 
    T1wBrainImage=${T1wFolder}/${T1wRestoreImageBrain}
    T1wBrainImageFile=`basename $T1wBrainImage`
    OutputTransform=${T1wFolder}/xfms/${fMRI2strOutputTransform} 
    BiasField=${T1wFolder}/${BiasField} 
    RegOutput=${fMRIFolder}/${RegOutput} 
    FreeSurferSubjectFolder=${i}/Preprocess_Structure/Freesurfer
    FreeSurferSubjectID=${Subject} 
    QAImage=${fMRIFolder}/${QAImage} 
    JacobianOut=${fMRIFolder}/${JacobianOut}

    if [ ! -f ${DCFolder}/${T1wBrainImageFile}.nii.gz ] ; then
        cp ${T1wBrainImage}.nii.gz ${DCFolder}/${T1wBrainImageFile}.nii.gz
    fi

    if [ ! -f ${QAImage}.nii.gz ] ; then
        echo "${i} BEGINNING EPI DISTORTION CORRECTION AND EPI TO T1W REGISTRATION"
    fi

    #### Use topup to distortion correct the scout scans (from TopupPreprocessingAll.sh, lines 165~412) ####
    # using a blip-reversed SE pair "fieldmap" sequence #
    topupdir=${DCFolder}/FieldMap
    if [ ! -d ${topupdir} ] ; then
        mkdir ${topupdir}
    fi
    PhaseEncodeOne=${SpinEchoPhaseEncodeNegative} 
    PhaseEncodeTwo=${SpinEchoPhaseEncodePositive} 
    DistortionCorrectionWarpFieldOutput=${DCFolder}/WarpField 
    DistortionCorrectionMagnitudeOutput=${DCFolder}/Magnitude    
    DistortionCorrectionMagnitudeBrainOutput=${DCFolder}/Magnitude_brain 
    DistortionCorrectionFieldOutput=${DCFolder}/TopupField  
    JacobianOutput=${DCFolder}/Jacobian 
    ## check dimensions of phase versus sbref images ##
    if [[ `fslhd ${PhaseEncodeOne} | grep '^dim[123]'` != `fslhd ${ScoutInputName} | grep '^dim[123]'` ]]
    then
        echo "${i} Resting Functional Readout Distortion Correction : ERROR Spin echo "fieldmap" has different dimensions than scout image, this requires a manual fix"
    fi
    ## check that the spin echo images match ##
    if [[ `fslhd ${PhaseEncodeOne} | grep '^dim[123]'` != `fslhd ${PhaseEncodeTwo} | grep '^dim[123]'` ]]
    then
        echo "${i} Resting Functional Readout Distortion Correction : ERROR Spin echo fieldmap images have different dimensions!"
    fi
    ## copy input images ##
    if [ ! -f ${topupdir}/PhaseOne.nii.gz ] || [ ! -f ${topupdir}/PhaseTwo.nii.gz ] || [ ! -f ${topupdir}/SBRef.nii.gz ] ; then
        ${FSLDIR}/bin/imcp ${PhaseEncodeOne} ${topupdir}/PhaseOne.nii.gz
        ${FSLDIR}/bin/imcp ${PhaseEncodeTwo} ${topupdir}/PhaseTwo.nii.gz
        ${FSLDIR}/bin/imcp ${ScoutInputName} ${topupdir}/SBRef.nii.gz
    fi

    ## Apply gradient non-linearity distortion correction to input images (SE pair) ##
    if [ ! -f ${topupdir}/PhaseOne_gdc.nii.gz ] || [ ! -f ${topupdir}/PhaseTwo_gdc.nii.gz ] || [ ! -f ${topupdir}/BothPhases.nii.gz ] ; then
        if [ ! $GradientDistortionCoeffs = "NONE" ] ; then
            echo "${i} Resting Functional Readout Distrotion Correction: Applying Gradient Non-linearity Distortion Correction to SE Pair"
            ${GlobalScripts}/GradientDistortionUnwarp.sh --workingdir=${topupdir} --coeffs=${GradientDistortionCoeffs} --in=${topupdir}/PhaseOne --out=${topupdir}/PhaseOne_gdc --owarp=${topupdir}/PhaseOne_gdc_warp
            ${GlobalScripts}/GradientDistortionUnwarp.sh --workingdir=${topupdir} --coeffs=${GradientDistortionCoeffs} --in=${topupdir}/PhaseTwo --out=${topupdir}/PhaseTwo_gdc --owarp=${topupdir}/PhaseTwo_gdc_warp
            if [[ $UseJacobian == "true" ]] ;then
                ${FSLDIR}/bin/fslmaths ${topupdir}/PhaseOne_gdc -mul ${topupdir}/PhaseOne_gdc_warp_jacobian ${topupdir}/PhaseOne_gdc
                ${FSLDIR}/bin/fslmaths ${topupdir}/PhaseTwo_gdc -mul ${topupdir}/PhaseTwo_gdc_warp_jacobian ${topupdir}/PhaseTwo_gdc
            fi      #overwrites inputs, no else needed
            #in the below stuff, the jacobians for both phases and sbref are applied unconditionally to a separate _jac image
            ## Make a dilated mask in the distortion corrected space ##
            ${FSLDIR}/bin/fslmaths ${topupdir}/PhaseOne -abs -bin -dilD ${topupdir}/PhaseOne_mask
            ${FSLDIR}/bin/applywarp --rel --interp=nn -i ${topupdir}/PhaseOne_mask -r ${topupdir}/PhaseOne_mask -w ${topupdir}/PhaseOne_gdc_warp -o ${topupdir}/PhaseOne_mask_gdc
            ${FSLDIR}/bin/fslmaths ${topupdir}/PhaseTwo -abs -bin -dilD ${topupdir}/PhaseTwo_mask
            ${FSLDIR}/bin/applywarp --rel --interp=nn -i ${topupdir}/PhaseTwo_mask -r ${topupdir}/PhaseTwo_mask -w ${topupdir}/PhaseTwo_gdc_warp -o ${topupdir}/PhaseTwo_mask_gdc
            ## Make a conservative (eroded) intersection of the two masks
            ${FSLDIR}/bin/fslmaths ${topupdir}/PhaseOne_mask_gdc -mas ${topupdir}/PhaseTwo_mask_gdc -ero -bin ${topupdir}/Mask
            ## Merge both sets of images ##
            ${FSLDIR}/bin/fslmerge -t ${topupdir}/BothPhases ${topupdir}/PhaseOne_gdc ${topupdir}/PhaseTwo_gdc
        else
            echo "${i} Resting Functional Readout Distrotion Correction: Not Applying Gradient Non-linearity Distortion Dorrection to SE Pair (but input names changed to match GDC files)"
            cp ${topupdir}/PhaseOne.nii.gz ${topupdir}/PhaseOne_gdc.nii.gz
            cp ${topupdir}/PhaseTwo.nii.gz ${topupdir}/PhaseTwo_gdc.nii.gz
            fslmerge -t ${topupdir}/BothPhases ${topupdir}/PhaseOne_gdc ${topupdir}/PhaseTwo_gdc
            fslmaths ${topupdir}/PhaseOne_gdc.nii.gz -mul 0 -add 1 ${topupdir}/Mask
        fi
    elif [ -e ${topupdir}/PhaseOne_gdc.nii.gz ] && [ -e ${topupdir}/PhaseTwo_gdc.nii.gz ] ; then
        echo "${i} Resting Functional Readout Distrotion Correction: Gradient Non-linearity Distortion Corrected on SE Pair (Or input names changed as if)"
    else
        echo "${i} Resting Functional Readout Distrotion Correction: ERROR Applying Gradient Non-linearity Distortion Dorrection to SE Pair (or changing input names to match GDC files)"
    fi

    ## Set up text files with all necessary parameters ##
    txtfname=${topupdir}/acqparams.txt
    dimtOne=`${FSLDIR}/bin/fslval ${topupdir}/PhaseOne dim4`
    dimtTwo=`${FSLDIR}/bin/fslval ${topupdir}/PhaseTwo dim4`
    if [ ! -f $txtfname ] ; then
        echo "${i} Resting Functional Readout Distortion Correction : Calculating the readout time and populating the parameter file appropriately"
        if [[ $UnwarpDir = "x" || $UnwarpDir = "x-" || $UnwarpDir = "-x" ]] ; then
            dimx=`${FSLDIR}/bin/fslval ${topupdir}/PhaseOne dim1`
            nPEsteps=$(($dimx - 1))
            ro_time=`echo "scale=6; ${DwellTime} * (${nPEsteps}/${grappa_factor})" | bc -l`     
            echo "      Total readout time is $ro_time secs"     
            j=1
            while [ $j -le $dimtOne ] ; do
                echo "-1 0 0 $ro_time" >> $txtfname
                ShiftOne="x-"
                j=`echo "$j + 1" | bc`
            done
            j=1
            while [ $j -le $dimtTwo ] ; do
                echo "1 0 0 $ro_time" >> $txtfname
                ShiftTwo="x"
                j=`echo "$j + 1" | bc`
            done
#       elif [[ $UnwarpDir = "y" || $UnwarpDir = "y-" || $UnwarpDir = "-y" ]] ; then
#           nPEsteps=$(($dimy - 1))
#           #Total_readout=DwellTime*(#of_PE_steps-1)
#           ro_time=`echo "scale=6; ${DwellTime} * ${nPEsteps}" | bc -l` #Compute Total_readout in secs with up to 6 decimal places
#           i=1
#           while [ $i -le $dimtOne ] ; do
#               echo "0 -1 0 $ro_time" >> $txtfname
#               ShiftOne="y-"
#               i=`echo "$i + 1" | bc`
#           done
#           i=1
#           while [ $i -le $dimtTwo ] ; do
#               echo "0 1 0 $ro_time" >> $txtfname
#               ShiftTwo="y"
#               i=`echo "$i + 1" | bc`
#           done
        fi
    else
        correct_lines=6
        check_lines=`cat ${txtfname} | wc -l`
        if [ $correct_lines -eq $check_lines ] ; then
            echo "${i} Resting Functional Readout Distortion Correction : Calculated the readout time and Created the parameter file appropriately"
        else
            echo "${i} Resting Functional Readout Distortion Correction : ERROR Calculating the readout time and Creating the parameter file appropriately"
        fi
    fi
    ## Pad in Z by one slice if odd so that topup does not complain (slice consists of zeros that will be dilated by following step) ##
    numslice=`fslval ${topupdir}/BothPhases dim3`
    if [ ! $(($numslice % 2)) -eq "0" ] ; then
        echo "      Padding Z by one slice"
        for Image in ${topupdir}/BothPhases ${topupdir}/Mask ; do
            fslroi ${Image} ${topupdir}/slice.nii.gz 0 -1 0 -1 0 1 0 -1
            fslmaths ${topupdir}/slice.nii.gz -mul 0 ${topupdir}/slice.nii.gz
            fslmerge -z ${Image} ${Image} ${topupdir}/slice.nii.gz
            rm ${topupdir}/slice.nii.gz
        done
    else
        echo "      Slice Number Even"
    fi

    ## Run Topup ##
    if [ ! -f ${topupdir}/WarpField_01.nii.gz ] ; then
        echo "${i} Resting Functional Readout Distortion Correction : Running Topup"
        ## Extrapolate the existing values beyond the mask (adding 1 just to avoid smoothing inside the mask) ##
        ${FSLDIR}/bin/fslmaths ${topupdir}/BothPhases -abs -add 1 -mas ${topupdir}/Mask -dilM -dilM -dilM -dilM -dilM ${topupdir}/BothPhases
        ## Topup ## 
        ${FSLDIR}/bin/topup --imain=${topupdir}/BothPhases --datain=$txtfname --config=${TopupConfig} --out=${topupdir}/Coefficents --iout=${topupdir}/Magnitudes --fout=${topupdir}/TopupField --dfout=${topupdir}/WarpField --rbmout=${topupdir}/MotionMatrix --jacout=${topupdir}/Jacobian -v
    elif [ -e ${topupdir}/WarpField_01.nii.gz ] && [ -e ${topupdir}/Coefficents_fieldcoef.nii.gz ] && [ -e ${topupdir}/Magnitudes.nii.gz ] && [ -e ${topupdir}/TopupField.nii.gz ] && [ -e ${topupdir}/MotionMatrix_01.mat ] && [ -e ${topupdir}/Jacobian_01.nii.gz ] ; then
        correct_topup=22
        coeffs=`ls ${topupdir}/Coefficents* | wc -l`
        jacs=`ls ${topupdir}/Jacobian_0* | wc -l`
        mag=`ls ${topupdir}/Magnitudes* | wc -l`
        motion_mats=`ls ${topupdir}/MotionMatrix_0* | wc -l`
        topupfield=`ls ${topupdir}/TopupField* | wc -l`
        warps=`ls ${topupdir}/WarpField_0* | wc -l`
        check_topup=$(( $coeffs + $jacs + $mag + $motion_mats + $topupfield + $warps )) 
        if [ $correct_topup -eq $check_topup ] ; then
            echo "${i} Resting Functional Readout Distortion Correction : Topup Completed"
        else
            echo "${i} Resting Functional Readout Distortion Correction : ERROR Running Topup"
        fi
    else
        echo "${i} Resting Functional Readout Distortion Correction : ERROR Running Topup"
    fi
    ## Remove Z slice padding if needed ##
    if [ ! $(($numslice % 2)) -eq "0" ] ; then
        echo "      Removing Z slice padding"
        for Image in ${topupdir}/BothPhases ${topupdir}/Mask ${topupdir}/Coefficents_fieldcoef ${topupdir}/Magnitudes ${topupdir}/TopupField* ${topupdir}/WarpField* ${topupdir}/Jacobian* ; do
            fslroi ${Image} ${Image} 0 -1 0 -1 0 ${numslice} 0 -1
        done
    fi

    ## Register SBRef (Scout) to Spin Echo EPI Image ##
    if [ ! -f ${topupdir}/WarpField.nii.gz ] || [ ! -f ${topupdir}/Jacobian.nii.gz ] ; then
        ## UNWARP DIR = x,y ##
        if [[ $UnwarpDir = "x" || $UnwarpDir = "y" ]] ; then
            echo "${i} Resting Functional Readout Distortion Correction : Registering Scout to PhaseTwo SE image (Unwarp Direction = $UnwarpDir)"
            ## select the first volume from PhaseTwo ##
            VolumeNumber=$(($dimtOne + 1))      #dimtOne=3
            vnum=`${FSLDIR}/bin/zeropad $VolumeNumber 2`        #vnum=04
            ## register scout to SE input (PhaseTwo) + combine motion and distortion correction ##
            ${FSLDIR}/bin/flirt -dof 6 -interp spline -in ${topupdir}/SBRef.nii.gz -ref ${topupdir}/PhaseTwo_gdc -omat ${topupdir}/SBRef2PhaseTwo_gdc.mat -out ${topupdir}/SBRef2PhaseTwo_gdc
            ${FSLDIR}/bin/convert_xfm -omat ${topupdir}/SBRef2WarpField.mat -concat ${topupdir}/MotionMatrix_${vnum}.mat ${topupdir}/SBRef2PhaseTwo_gdc.mat
            ${FSLDIR}/bin/convertwarp --relout --rel -r ${topupdir}/PhaseTwo_gdc --premat=${topupdir}/SBRef2WarpField.mat --warp1=${topupdir}/WarpField_${vnum} --out=${topupdir}/WarpField.nii.gz
            ${FSLDIR}/bin/imcp ${topupdir}/Jacobian_${vnum}.nii.gz ${topupdir}/Jacobian.nii.gz
            SBRefPhase=Two
        ## UNWARP DIR = -x,-y ##
        elif [[ $UnwarpDir = "x-" || $UnwarpDir = "-x" || $UnwarpDir = "y-" || $UnwarpDir = "-y" ]] ; then
            echo "${i} Resting Functional Readout Distortion Correction : Registering Scout to PhaseOne SE image (Unwarp Direction = $UnwarpDir)"
            ## select the first volume from PhaseOne ##
            VolumeNumber=$((0 + 1))
            vnum=`${FSLDIR}/bin/zeropad $VolumeNumber 2`
            ## register scout to SE input (PhaseOne) + combine motion and distortion correction ##
            ${FSLDIR}/bin/flirt -dof 6 -interp spline -in ${topupdir}/SBRef.nii.gz -ref ${topupdir}/PhaseOne_gdc -omat ${topupdir}/SBRef2PhaseOne_gdc.mat -out ${topupdir}/SBRef2PhaseOne_gdc
            ${FSLDIR}/bin/convert_xfm -omat ${topupdir}/SBRef2WarpField.mat -concat ${topupdir}/MotionMatrix_${vnum}.mat ${topupdir}/SBRef2PhaseOne_gdc.mat
            ${FSLDIR}/bin/convertwarp --relout --rel -r ${topupdir}/PhaseOne_gdc --premat=${topupdir}/SBRef2WarpField.mat --warp1=${topupdir}/WarpField_${vnum} --out=${topupdir}/WarpField.nii.gz
            ${FSLDIR}/bin/imcp ${topupdir}/Jacobian_${vnum}.nii.gz ${topupdir}/Jacobian.nii.gz
            SBRefPhase=One
        fi
    elif [ -e ${topupdir}/WarpField.nii.gz ] && [ -e ${topupdir}/Jacobian.nii.gz ] ; then
        echo "${i} Resting Functional Readout Distortion Correction : Scout Registered to SE image (Unwarp Direction = $UnwarpDir)"
    else
        echo "${i} Resting Functional Readout Distortion Correction : ERROR Registering Scout to SE image (Unwarp Direction = $UnwarpDir)"
    fi

    ## Distortion Correction ##
    if [ ! -f ${topupdir}/PhaseTwo_gdc_dc.nii.gz ] || [ ! -f ${topupdir}/PhaseTwo_gdc_dc_jac.nii.gz ] ; then
        echo "${i} Resting Functional Readout Distortion Correction : Correcting PhaseTwo (BLIP_RL) image" 
        # PhaseTwo (first vol) - warp and Jacobian modulate to get distortion corrected output
        VolumeNumber=$(($dimtOne + 1))
        vnum=`${FSLDIR}/bin/zeropad $VolumeNumber 2`
        ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${topupdir}/PhaseTwo_gdc -r ${topupdir}/PhaseTwo_gdc --premat=${topupdir}/MotionMatrix_${vnum}.mat -w ${topupdir}/WarpField_${vnum} -o ${topupdir}/PhaseTwo_gdc_dc
        ${FSLDIR}/bin/fslmaths ${topupdir}/PhaseTwo_gdc_dc -mul ${topupdir}/Jacobian_${vnum} ${topupdir}/PhaseTwo_gdc_dc_jac
    elif [ -e ${topupdir}/PhaseTwo_gdc_dc.nii.gz ] && [ -e ${topupdir}/PhaseTwo_gdc_dc_jac.nii.gz ] ; then
        echo "${i} Resting Functional Readout Distortion Correction : PhaseTwo (BLIP_RL) Image Corrected" 
    else
        echo "${i} Resting Functional Readout Distortion Correction : ERROR Correcting PhaseTwo (BLIP_RL) image" 
    fi
    if [ ! -f ${topupdir}/PhaseOne_gdc_dc.nii.gz ] || [ ! -f ${topupdir}/PhaseOne_gdc_dc_jac.nii.gz ] ; then
        echo "${i} Resting Functional Readout Distortion Correction : Correcting PhaseOne (BLIP_LR) image" 
        # PhaseOne (first vol) - warp and Jacobian modulate to get distortion corrected output
        VolumeNumber=$((0 + 1))
        vnum=`${FSLDIR}/bin/zeropad $VolumeNumber 2`
        ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${topupdir}/PhaseOne_gdc -r ${topupdir}/PhaseOne_gdc --premat=${topupdir}/MotionMatrix_${vnum}.mat -w ${topupdir}/WarpField_${vnum} -o ${topupdir}/PhaseOne_gdc_dc
        ${FSLDIR}/bin/fslmaths ${topupdir}/PhaseOne_gdc_dc -mul ${topupdir}/Jacobian_${vnum} ${topupdir}/PhaseOne_gdc_dc_jac
    elif [ -e ${topupdir}/PhaseOne_gdc_dc.nii.gz ] && [ -e ${topupdir}/PhaseOne_gdc_dc_jac.nii.gz ] ; then
        echo "${i} Resting Functional Readout Distortion Correction : PhaseOne (BLIP_LR) Image Corrected" 
    else
        echo "${i} Resting Functional Readout Distortion Correction : ERROR Correcting PhaseOne (BLIP_LR) image" 
    fi
    if [ ! -f ${topupdir}/SBRef_dc.nii.gz ] || [ ! -f ${topupdir}/SBRef_dc_jac.nii.gz ] ; then
        echo "${i} Resting Functional Readout Distortion Correction : Correcting Scout (SBRef) image" 
        # Scout - warp and Jacobian modulate to get distortion corrected output
        ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${topupdir}/SBRef.nii.gz -r ${topupdir}/SBRef.nii.gz -w ${topupdir}/WarpField.nii.gz -o ${topupdir}/SBRef_dc.nii.gz
        ${FSLDIR}/bin/fslmaths ${topupdir}/SBRef_dc.nii.gz -mul ${topupdir}/Jacobian.nii.gz ${topupdir}/SBRef_dc_jac.nii.gz
    elif [ -e ${topupdir}/SBRef_dc.nii.gz ] && [ -e ${topupdir}/SBRef_dc_jac.nii.gz ] ; then
        echo "${i} Resting Functional Readout Distortion Correction : Scout (SBRef) Image Corrected"
    else
        echo "${i} Resting Functional Readout Distortion Correction : ERROR Correcting Scout (SBRef) image" 
    fi

    ## Calculate Equivalent Field Map ##
    if [ ! -f ${topupdir}/Magnitude_brain.nii.gz ] ; then
        echo "${i} Resting Functional Readout Distortion Correction : Calculating Equivalent Field Map"
        ${FSLDIR}/bin/fslmaths ${topupdir}/TopupField -mul 6.283 ${topupdir}/TopupField
        ${FSLDIR}/bin/fslmaths ${topupdir}/Magnitudes.nii.gz -Tmean ${topupdir}/Magnitude.nii.gz
        ## Brain extract the magnitude image ##
        ${FSLDIR}/bin/bet ${topupdir}/Magnitude ${topupdir}/Magnitude_brain -f 0.35 -m
    else
        echo "${i} Resting Functional Readout Distortion Correction : Equivalent Field Map Calculated"
    fi

    ## copy images to specified outputs ##
    if [ ! -f ${DistortionCorrectionWarpFieldOutput}.nii.gz ] ; then
        ${FSLDIR}/bin/imcp ${topupdir}/WarpField.nii.gz ${DistortionCorrectionWarpFieldOutput}.nii.gz
        ${FSLDIR}/bin/imcp ${topupdir}/Jacobian.nii.gz ${JacobianOutput}.nii.gz
        ${FSLDIR}/bin/imcp ${topupdir}/TopupField.nii.gz ${DistortionCorrectionFieldOutput}.nii.gz
        ${FSLDIR}/bin/imcp ${topupdir}/Magnitude.nii.gz ${DistortionCorrectionMagnitudeOutput}.nii.gz
        ${FSLDIR}/bin/imcp ${topupdir}/Magnitude_brain.nii.gz ${DistortionCorrectionMagnitudeBrainOutput}.nii.gz
    elif [ -e ${DistortionCorrectionWarpFieldOutput}.nii.gz ] && [ -e ${JacobianOutput}.nii.gz ] && [ -e ${DistortionCorrectionFieldOutput}.nii.gz ] && [ -e ${DistortionCorrectionMagnitudeOutput}.nii.gz ] && [ -e ${DistortionCorrectionMagnitudeBrainOutput}.nii.gz ] ; then
        echo "${i} Resting Functional Readout Distortion Correction Completed"
    else
        echo "${i} Resting Functional Readout Distortion Correction Not Completed"
    fi
    #### COMPLETED TOPUP RELATED PROCESSES (lines 165~412) ####

    #### Create a spline interpolated image of scout (distortion corrected in same space) ####
    if [ ! -f ${DCFolder}/${ScoutInputFile}_undistorted.nii.gz ] ; then
        echo "${i} Resting Functional Pre-Register : Creating a spline interpolated image of scout (distortion corrected in same space)"
        ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${ScoutInputName} -r ${ScoutInputName} -w ${DCFolder}/WarpField.nii.gz -o ${DCFolder}/${ScoutInputFile}_undistorted
        # gdc jacobian is already applied in main script, where the gdc call for the scout is
        if [[ $UseJacobian == "true" ]] ; then
            echo "      Applying Jacobian correction to scout image"
            ${FSLDIR}/bin/fslmaths ${DCFolder}/${ScoutInputFile}_undistorted -mul ${DCFolder}/Jacobian.nii.gz ${DCFolder}/${ScoutInputFile}_undistorted
        fi
    else
        echo "${i} Resting Functional Pre-Register : Spline Interpolated Image of Scout Created (distortion corrected in same space)"
    fi

    #### Register undistorted scout image to T1w (from epi_reg_dof) ####
    # this is just an initial registration, refined later in this script, but it is actually pretty good
    vepi=${DCFolder}/${ScoutInputFile}_undistorted
    vrefhead=${T1wImage}
    vrefbrain=${DCFolder}/${T1wBrainImageFile}
    vout=${DCFolder}/${ScoutInputFile}_undistorted2T1w_init
    if [ ! -f ${vout}.nii.gz ]
    then
        echo "${i} Resting Functional Initial Registration : Register undistorted scout image to T1w (epi_reg_dof)"    
        ## create the WM segmentation ##
        if [ ! -f ${vout}_fast_wmseg ]
        then
            echo "      Running FAST segmentation"
            $FSLDIR/bin/fast -o ${vout}_fast ${vrefbrain}
            $FSLDIR/bin/fslmaths ${vout}_fast_pve_2 -thr 0.5 -bin ${vout}_fast_wmseg
        else
            echo "      FAST segmentation completed"
        fi
        ## make a WM edge map for visualisation (good to overlay in FSLView) ##
        if [ ! -f ${vout}_fast_wmedge ]
        then
            $FSLDIR/bin/fslmaths ${vout}_fast_wmseg -edge -bin -mas ${vout}_fast_wmseg ${vout}_fast_wmedge
        else
            echo "      WM edge for visualisation created"
        fi
        ## do a standard flirt pre-alignment ##
        if [ ! -f ${vout}_init.mat ]
        then
            echo "      Running FLIRT pre-alignment"
            $FSLDIR/bin/flirt -ref ${vrefbrain} -in ${vepi} -dof ${dof} -omat ${vout}_init.mat
        else
            echo "      FLIRT pre-alignment Completed"
        fi
        ## Run BBR ##
        echo "      Running BBR"
        $FSLDIR/bin/flirt -ref ${vrefhead} -in ${vepi} -dof ${dof} -cost bbr -wmseg ${vout}_fast_wmseg -init ${vout}_init.mat -omat ${vout}.mat -out ${vout} -schedule ${FSLDIR}/etc/flirtsch/bbr.sch
        $FSLDIR/bin/applywarp -i ${vepi} -r ${vrefhead} -o ${vout} --premat=${vout}.mat --interp=spline
        cp ${vout}.mat ${DCFolder}/fMRI2str.mat
    else
        echo "${i} Resting Functional Initial Registration : Undistorted Scout Image Registered to T1w (epi_reg_dof Completed)"
    fi
    #### COMPLETED epi_reg_dof (lines 429~469) ####

    if [ ! -f ${DCFolder}/fMRI2str.mat ] ; then
        cp ${DCFolder}/${ScoutInputFile}_undistorted2T1w_init.mat ${DCFolder}/fMRI2str.mat
    fi

    #### Generate combined warpfields and spline interpolated images + apply bias field correction ####
    if [ ! -f ${DCFolder}/${ScoutInputFile}_undistorted2T1w_init_warp.nii.gz ] ; then
        echo "${i} Resting Functional Initial Registration : Generating combined warpfields and spline interpolated images and Apply bias field correction"
        ${FSLDIR}/bin/convertwarp --relout --rel -r ${T1wImage} --warp1=${DCFolder}/WarpField.nii.gz --postmat=${DCFolder}/${ScoutInputFile}_undistorted2T1w_init.mat -o ${DCFolder}/${ScoutInputFile}_undistorted2T1w_init_warp
        ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${DCFolder}/Jacobian.nii.gz -r ${T1wImage} --premat=${DCFolder}/${ScoutInputFile}_undistorted2T1w_init.mat -o ${DCFolder}/Jacobian2T1w.nii.gz
    elif [ -e ${DCFolder}/${ScoutInputFile}_undistorted2T1w_init_warp.nii.gz ] && [ -e ${DCFolder}/Jacobian2T1w.nii.gz ] ; then
        echo "${i} Resting Functional Initial Registration : Combined Warpfields and Spline Interpolated Images Generated and Bias Field Correction Applied"
    else
        echo "${i} Resting Functional Initial Registration : ERROR Generating combined warpfields and spline interpolated images and Apply bias field correction"
    fi

    #### 1-Step Resample from input (gdc) scout - NOTE: no longer includes jacobian correction, if specified ####
    if [ ! -f ${DCFolder}/${ScoutInputFile}_undistorted2T1w_init_init.nii.gz ] ; then
        echo "${i} Resting Functional Initial Registration : 1-step Resampling from input scout"
        mv ${DCFolder}/${ScoutInputFile}_undistorted2T1w_init.nii.gz ${DCFolder}/${ScoutInputFile}_undistorted2T1w_init_init.nii.gz
        ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${ScoutInputName} -r ${T1wImage} -w ${DCFolder}/${ScoutInputFile}_undistorted2T1w_init_warp -o ${DCFolder}/${ScoutInputFile}_undistorted2T1w_init
    else
        echo "${i} Resting Functional Initial Registration : 1-Step Resampling from Input Scout Completed"
    fi

    ##### Resample phase images to T1w space ####
    Files="PhaseOne_gdc_dc PhaseTwo_gdc_dc SBRef_dc"
    for File in ${Files} ; do
        if [ ! -f ${DCFolder}/$File.nii.gz ] ; then
            echo "${i} Resting Functional Initial Registration : Resampling $File to T1w space"
            if [[ $UseJacobian == "true" ]] ; then
                ${FSLDIR}/bin/applywarp --interp=spline -i ${DCFolder}/FieldMap/${File}_jac -r ${T1wFolder}/T2w_acpc_dc.nii.gz --premat=${DCFolder}/fMRI2str.mat -o ${DCFolder}/${File}
            else
                ${FSLDIR}/bin/applywarp --interp=spline -i ${DCFolder}/FieldMap/${File} -r ${T1wFolder}/T2w_acpc_dc.nii.gz --premat=${DCFolder}/fMRI2str.mat -o ${DCFolder}/${File}
            fi
        elif [ -e ${DCFolder}/$File.nii.gz ] ; then
            echo "${i} Resting Functional Initial Registration : Resampled $File to T1w space"
        else
            echo "${i} Resting Functional Initial Registration : ERROR Resampling $File to T1w space"
        fi
    done

    #### Bias Correction ####
    if [ "$BiasCorrection" == "SEBASED" ] ; then
        echo "${i} Resting Functional Initial Registration : Performing $BiasCorrection bias correction"
        mkdir -p $DCFolder/ComputeSpinEchoBiasField
        bash ${HCPPIPEDIR_fMRIVol}/edit_ComputeSpinEchoBiasField.sh --workingdir=$DCFolder/ComputeSpinEchoBiasField --subjectfolder=$SubjectFolder --fmriname=$NamingOffMRI --corticallut=$HCPPIPEDIR/global/config/FreeSurferCorticalLabelTableLut.txt --subcorticallut=$HCPPIPEDIR/global/config/FreeSurferSubcorticalLabelTableLut.txt --smoothingfwhm=2 --inputdir=$DCFolder
    fi

    case "$BiasCorrection" in
    NONE)
        UseBiasField=""
        UseBiasFieldMNI=""
        ;;
    LEGACY)
        UseBiasField="${BiasField}"
        UseBiasFieldMNI="${fMRIFolder}/${BiasFieldMNI}.${FinalfMRIResolution}"
        ;;
    SEBASED)
        UseBiasField="${DCFolder}/ComputeSpinEchoBiasField/${NamingOffMRI}_sebased_bias.nii.gz"
        UseBiasFieldMNI="$sebasedBiasFieldMNI"
        ;;
    esac

    #### Apply Jacobian correction and bias correction options to scout image ####
    if [ ! -f ${DCFolder}/${ScoutInputFile}_undistorted2T1w_init_preBC.nii.gz ] ; then
        if [[ $UseJacobian == "true" ]] ; then
            echo "${i} Resting Functional Initial Registration : Applying Jacobian correction and bias correction options to scout image"
            cp ${DCFolder}/${ScoutInputFile}_undistorted2T1w_init.nii.gz ${DCFolder}/${ScoutInputFile}_undistorted2T1w_init_preBC.nii.gz
            if [[ "$UseBiasField" != "" ]] ; then
                ${FSLDIR}/bin/fslmaths ${DCFolder}/${ScoutInputFile}_undistorted2T1w_init_preBC -div ${UseBiasField} -mul ${DCFolder}/Jacobian2T1w.nii.gz ${DCFolder}/${ScoutInputFile}_undistorted2T1w_init.nii.gz
            else
                ${FSLDIR}/bin/fslmaths ${DCFolder}/${ScoutInputFile}_undistorted2T1w_init_preBC -mul ${DCFolder}/Jacobian2T1w.nii.gz ${DCFolder}/${ScoutInputFile}_undistorted2T1w_init.nii.gz
            fi
        else
            echo "${i} Resting Functional Initial Registration : Not Applying Jacobian correction and bias correction options to scout image"
            cp ${DCFolder}/${ScoutInputFile}_undistorted2T1w_init.nii.gz ${DCFolder}/${ScoutInputFile}_undistorted2T1w_init_preBC.nii.gz
            if [[ "$UseBiasField" != "" ]] ; then
                ${FSLDIR}/bin/fslmaths ${DCFolder}/${ScoutInputFile}_undistorted2T1w_init_preBC -div ${UseBiasField} ${DCFolder}/${ScoutInputFile}_undistorted2T1w_init.nii.gz
            fi
        fi
    elif [ -e ${DCFolder}/${ScoutInputFile}_undistorted2T1w_init_preBC.nii.gz ] && [ -e ${DCFolder}/${ScoutInputFile}_undistorted2T1w_init.nii.gz ] ; then
        echo "${i} Resting Functional Initial Registration : Applied Relevant Jacobian Correction and Bias Correction Options to Scout Image"
    else
        echo "${i} Resting Functional Initial Registration : ERROR Applying Relevant Jacobian Correction and Bias Correction Options to Scout Image"
    fi
    # Rplaced "1-step resampled" output (which used output from epi_reg_dof) 

    ### FREESURFER BBR - found to be an improvement, probably due to better GM/WM boundary ####
    if [ ! -f ${DCFolder}/${ScoutInputFile}_undistorted2T1w.nii.gz ] ; then
        echo "${i} Resting Functional Fine Tuning Registration : Running bbregister"
        SUBJECTS_DIR=${FreeSurferSubjectFolder}
        export SUBJECTS_DIR 
        ${FREESURFER_HOME}/bin/bbregister --s ${FreeSurferSubjectID} --mov ${DCFolder}/${ScoutInputFile}_undistorted2T1w_init.nii.gz --surf white.deformed --init-reg ${FreeSurferSubjectFolder}/${FreeSurferSubjectID}/mri/transforms/eye.dat --bold --reg ${DCFolder}/EPItoT1w.dat --${dof} --o ${DCFolder}/${ScoutInputFile}_undistorted2T1w.nii.gz
        ## Create FSL-style matrix and then combine with existing warp fields                                               
        if [ ! -f ${DCFolder}/fMRI2str_refinement.mat ] ; then
            echo "${i} Resting Functional Fine Tuning Registration : Creating FSL-style matrix and then combining with existing warp fields"
            ${FREESURFER_HOME}/bin/tkregister2 --noedit --reg ${DCFolder}/EPItoT1w.dat --mov ${DCFolder}/${ScoutInputFile}_undistorted2T1w_init.nii.gz --targ ${T1wImage}.nii.gz --fslregout ${DCFolder}/fMRI2str_refinement.mat
        fi
    elif [ -e ${DCFolder}/${ScoutInputFile}_undistorted2T1w.nii.gz ] && [ -e ${DCFolder}/fMRI2str_refinement.mat ] ; then
        echo "${i} Resting Functional Fine Tuning Registration : Completed bbregister, Created FSL-style Matrix and Combined with Existing Warp Fields "
    else
        echo "${i} Resting Functional Fine Tuning Registration : ERROR running bbregister, creatihg FSL-style matrix and combining with Existing Warp Fields "
    fi

    if [ ! -f ${DCFolder}/fMRI2str.nii.gz ] ; then
        ${FSLDIR}/bin/convertwarp --relout --rel --warp1=${DCFolder}/${ScoutInputFile}_undistorted2T1w_init_warp.nii.gz --ref=${T1wImage} --postmat=${DCFolder}/fMRI2str_refinement.mat --out=${DCFolder}/fMRI2str.nii.gz
    fi

    #### Create final affine from undistorted fMRI space to T1w space (will need it if it making SEBASED bias field) ####
    if [ ! -f ${DCFolder}/fMRI2str_init.mat ] ; then
        echo "${i} Resting Functional Fine Tuning Registration : Creating final affine from undistorted fMRI space to T1w space"
        cp ${DCFolder}/fMRI2str.mat ${DCFolder}/fMRI2str_init.mat       ## I added to keep file 
        ${FSLDIR}/bin/convert_xfm -omat ${DCFolder}/fMRI2str.mat -concat ${DCFolder}/fMRI2str_refinement.mat ${DCFolder}/${ScoutInputFile}_undistorted2T1w_init.mat
    elif [ -e ${DCFolder}/fMRI2str_init.mat ] && [ -e ${DCFolder}/fMRI2str.mat ] ; then
        echo "${i} Resting Functional Fine Tuning Registration : Created Final Affine from Undistorted fMRI Space to T1w Space"
    else
        echo "${i} Resting Functional Fine Tuning Registration : ERROR Creating final affine from undistorted fMRI space to T1w space"
    fi

    #### Resample SE field maps ####
    Files="PhaseOne_gdc_dc PhaseTwo_gdc_dc SBRef_dc"
    for File in ${Files} ; do
        if [ ! -f ${DCFolder}/${File}_init_resample.nii.gz ] ; then 
            echo "${i} Resting Functional Fine Tuning Registration : Resampling ${File}"
            ${FSLDIR}/bin/imcp ${DCFolder}/${File} ${DCFolder}/${File}_init_resample        ## I added to keep file 
            if [[ $UseJacobian == "true" ]] ; then
                ${FSLDIR}/bin/applywarp --interp=spline -i "${DCFolder}/FieldMap/${File}_jac" -r ${T1wFolder}/T2w_acpc_dc.nii.gz --premat=${DCFolder}/fMRI2str.mat -o ${DCFolder}/${File}
            else                           
                ${FSLDIR}/bin/applywarp --interp=spline -i "${DCFolder}/FieldMap/${File}" -r ${T1wFolder}/T2w_acpc_dc.nii.gz --premat=${DCFolder}/fMRI2str.mat -o ${DCFolder}/${File}
            fi
        elif [ -e ${DCFolder}/${File}_init_resample.nii.gz ] && [ -e ${DCFolder}/${File}.nii.gz ] ; then 
            echo "${i} Resting Functional Fine Tuning Registration : Resampled ${File}"
        else
            echo "${i} Resting Functional Fine Tuning Registration : ERROR Resampling ${File}"
        fi
    done
        
    #### Final bias field computation ####
    if [ ! -f ${DCFolder}/PhaseOne_gdc_dc_unbias.nii.gz ] || [ ! -f ${DCFolder}/PhaseTwo_gdc_dc_unbias.nii.gz ] ; then
        echo "${i} Resting Functional Fine Tuning Registration : Running Final bias field computation"
        if [[ $BiasCorrection == "SEBASED" ]] ; then
            # reusing the same working dir as previous run
            "${HCPPIPEDIR_fMRIVol}/ComputeSpinEchoBiasField.sh" --workingdir="$DCFolder/ComputeSpinEchoBiasField" --subjectfolder="$SubjectFolder" --fmriname="$NamingOffMRI" --corticallut="$HCPPIPEDIR/global/config/FreeSurferCorticalLabelTableLut.txt" --subcorticallut="$HCPPIPEDIR/global/config/FreeSurferSubcorticalLabelTableLut.txt" --smoothingfwhm="2" --inputdir="$DCFolder"                       
            # don't need to do anything more with scout, it is 1-step resampled and bias correction, jacobians reapplied
            Files="PhaseOne_gdc_dc PhaseTwo_gdc_dc"    
            for File in ${Files} ; do
                ## apply the new bias field to them for output
                ${FSLDIR}/bin/fslmaths ${DCFolder}/${File} -div "$UseBiasField" ${DCFolder}/${File}_unbias
            done
            ## copy bias field and dropouts, etc to results dir 
            ${FSLDIR}/bin/imcp "$DCFolder/ComputeSpinEchoBiasField/${NamingOffMRI}_dropouts" "$SubjectFolder/T1w/Results/$NamingOffMRI/${NamingOffMRI}_dropouts"
            ${FSLDIR}/bin/imcp "$DCFolder/ComputeSpinEchoBiasField/${NamingOffMRI}_sebased_bias" "$SubjectFolder/T1w/Results/$NamingOffMRI/${NamingOffMRI}_sebased_bias"
            ${FSLDIR}/bin/imcp "$DCFolder/ComputeSpinEchoBiasField/${NamingOffMRI}_sebased_reference" "$SubjectFolder/T1w/Results/$NamingOffMRI/${NamingOffMRI}_sebased_reference"
        else
        # don't need to do anything more with scout, it is 1-step resampled and bias correction, jacobians reapplied
            Files="PhaseOne_gdc_dc PhaseTwo_gdc_dc"
            for File in ${Files} ; do
                if [[ $UseBiasField ]] ; then
                    ## apply the bias field to them for output (really only the phase images, but whatever)
                    ${FSLDIR}/bin/fslmaths ${DCFolder}/${File} -div "$UseBiasField" ${DCFolder}/${File}_unbias
                else
                    ${FSLDIR}/bin/imcp ${DCFolder}/${File} ${DCFolder}/${File}_unbias
                fi
            done
        fi
    elif [ -e ${DCFolder}/PhaseOne_gdc_dc_unbias.nii.gz ] && [ -e ${DCFolder}/PhaseTwo_gdc_dc.nii.gz ] ; then
        echo "${i} Resting Functional Fine Tuning Registration : Final Bias Field Computation Completed"
    else
        echo "${i} Resting Functional Fine Tuning Registration : ERROR Running Final bias field computation"
    fi

    #### Create warped image with spline interpolation, bias correction and (optional) Jacobian modulation ####
    #NOTE: Jacobian2T1w should be only the topup or fieldmap warpfield's jacobian, not including the gdc warp
    #the input scout is the gdc scout, which should already have had the gdc jacobian applied by the main script
    if [ ! -f ${DCFolder}/${ScoutInputFile}_undistorted2T1w_bbregister.nii.gz ] ; then
        ${FSLDIR}/bin/imcp ${DCFolder}/${ScoutInputFile}_undistorted2T1w ${DCFolder}/${ScoutInputFile}_undistorted2T1w_bbregister       ## I added to keep file
        echo "${i} Resting Functional Fine Tuning Registration : Creating warped image with spline interpolation, bias correction and (optional) Jacobian modulation"
        ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${ScoutInputName} -r ${T1wImage}.nii.gz -w ${DCFolder}/fMRI2str.nii.gz -o ${DCFolder}/${ScoutInputFile}_undistorted2T1w
        ${FSLDIR}/bin/imcp ${DCFolder}/${ScoutInputFile}_undistorted2T1w ${DCFolder}/${ScoutInputFile}_undistorted2T1w_preRegOutput       ## I added to keep file
    elif [ -e ${DCFolder}/${ScoutInputFile}_undistorted2T1w_bbregister.nii.gz ] && [ -e ${DCFolder}/${ScoutInputFile}_undistorted2T1w.nii.gz ]  ; then
        echo "${i} Resting Functional Fine Tuning Registration : Created Warped Image with Spline Interpolation, Bias Correction and (optional) Jacobian Modulation"
    else
        echo "${i} Resting Functional Fine Tuning Registration : ERROR Creating warped image with spline interpolation, bias correction and (optional) Jacobian modulation"
    fi

    #### Resample fieldmap jacobian with new registration ####
    if [ ! -f ${DCFolder}/Jacobian2T1w_init.nii.gz ] ; then
        ${FSLDIR}/bin/imcp ${DCFolder}/Jacobian2T1w ${DCFolder}/Jacobian2T1w_init       ## I added to keep file
        echo "${i} Resting Functional Fine Tuning Registration : Resampling fieldmap jacobian with new registration"
        ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${DCFolder}/Jacobian.nii.gz -r ${T1wImage} --premat=${DCFolder}/fMRI2str.mat -o ${DCFolder}/Jacobian2T1w.nii.gz
        if [[ $UseJacobian == "true" ]] ; then
            echo "      applying Jacobian modulation"
            if [[ "$UseBiasField" != "" ]] ; then
                ${FSLDIR}/bin/fslmaths ${DCFolder}/${ScoutInputFile}_undistorted2T1w -div ${UseBiasField} -mul ${DCFolder}/Jacobian2T1w.nii.gz ${DCFolder}/${ScoutInputFile}_undistorted2T1w
            else
                ${FSLDIR}/bin/fslmaths ${DCFolder}/${ScoutInputFile}_undistorted2T1w -mul ${DCFolder}/Jacobian2T1w.nii.gz ${DCFolder}/${ScoutInputFile}_undistorted2T1w
            fi
        else
            echo "      Not applying Jacobian modulation"
            if [[ "$UseBiasField" != "" ]] ; then
                echo "      Using Bias Field ${UseBiasField}"
                ${FSLDIR}/bin/fslmaths ${DCFolder}/${ScoutInputFile}_undistorted2T1w -div ${UseBiasField} ${DCFolder}/${ScoutInputFile}_undistorted2T1w
            fi
        fi
    elif [ -e ${DCFolder}/Jacobian2T1w_init.nii.gz ] && [ -e ${DCFolder}/Jacobian2T1w.nii.gz ] && [ -e ${DCFolder}/${ScoutInputFile}_undistorted2T1w_preRegOutput.nii.gz ] && [ -e ${DCFolder}/${ScoutInputFile}_undistorted2T1w.nii.gz ] ; then
        echo "${i} Resting Functional Fine Tuning Registration : Resampled Fieldmap Jacobian with New Registration"
    else
        echo "${i} Resting Functional Fine Tuning Registration : ERROR Resampling fieldmap jacobian with new registration"
    fi

    #### Sorting outputs of EPI Distortion Correction and EPI to T1w Registration ####
    if [ ! -f ${RegOutput}.nii.gz ] ; then
        cp ${DCFolder}/${ScoutInputFile}_undistorted2T1w.nii.gz ${RegOutput}.nii.gz
    fi
    OutputTransformDir=$(dirname ${OutputTransform})
    if [ ! -e ${OutputTransformDir} ] ; then
        mkdir -p ${OutputTransformDir}
    fi
    if [ ! -f ${OutputTransform}.nii.gz ] ; then
        cp ${DCFolder}/fMRI2str.nii.gz ${OutputTransform}.nii.gz
    fi
    if [ ! -f ${JacobianOut}.nii.gz ] ; then
        cp ${DCFolder}/Jacobian2T1w.nii.gz ${JacobianOut}.nii.gz
    fi
    ## QA image (sqrt of EPI * T1w)
    if [ ! -f ${QAImage}.nii.gz ] ; then
        ${FSLDIR}/bin/fslmaths ${T1wRestoreImage}.nii.gz -mul ${RegOutput}.nii.gz -sqrt ${QAImage}.nii.gz
    fi

    if [ -e ${RegOutput}.nii.gz ] && [ -e ${OutputTransform}.nii.gz ] && [ -e ${JacobianOut}.nii.gz ] && [ -e ${QAImage}.nii.gz ] ; then
        echo "${i} COMPLETED EPI DISTORTION CORRECTION AND EPI TO T1W REGISTRATION"     # lines 136~699
    fi

#### One Step Resampling (from OneStepResampling.sh) ####
    onestep=${fMRIFolder}/OneStepResampling
    if [ ! -d ${onestep} ] ; then
        mkdir ${onestep}
    fi
    InputfMRI2=${fMRIFolder}/${OrigTCSName}.nii.gz        
    T1wImage2=${AtlasSpaceFolder}/${T1wAtlasName}        
    StructuralToStandard=${AtlasSpaceFolder}/xfms/${AtlasTransform}
    OutputTransform2=${AtlasSpaceFolder}/xfms/${OutputfMRI2StandardTransform}
    OutputInvTransform=${AtlasSpaceFolder}/xfms/${Standard2OutputfMRITransform}
    MotionMatrixFolder_mcflirt=${fMRIFolder}/${MotionMatrixFolder}_mcflirt 
    MotionMatrixFolder_flirt=${fMRIFolder}/${MotionMatrixFolder}_flirt 
    OutputfMRI=${fMRIFolder}/${NamingOffMRI}_nonlin 
    FreeSurferBrainMask=${AtlasSpaceFolder}/${FreeSurferBrainMask}
    #BiasField2=${AtlasSpaceFolder}/${BiasFieldMNI}
    BiasField2=${AtlasSpaceFolder}/${BiasFieldMNI}      ## REQUIRES POSTFREESURFER
    GradientDistortionField=${fMRIFolder}/${NamingOffMRI}_gdc_warp
    ScoutInput=${fMRIFolder}/${OrigScoutName}   
    ScoutOutput=${fMRIFolder}/${NamingOffMRI}_SBRef_nonlin 
    JacobianOut2=${fMRIFolder}/Jacobian_MNI.${FinalfMRIResolution} 
    BiasFieldFile2=`basename "$BiasField2"`
    T1wImageFile2=`basename $T1wImage2`
    FreeSurferBrainMaskFile=`basename "$FreeSurferBrainMask"`

    if [ ! -e ${ScoutOutput}.nii.gz ] ; then
        echo "${i} RESTING FUNCTIONAL ONE STEP RESAMPLING"
        ## Save TR for later
        TR_vol=`${FSLDIR}/bin/fslval ${InputfMRI2} pixdim4 | cut -d " " -f 1`
        echo "TR=$TR_vol"
        NumFrames=`${FSLDIR}/bin/fslval ${InputfMRI2} dim4`
        echo "NumFrames=$NumFrames"

        #### Create fMRI resolution standard space files for T1w image, wmparc, and brain mask ####
        #NB: don't use FLIRT to do spline interpolation with -applyisoxfm for the 2mm and 1mm cases because it doesn't know the peculiarities of the MNI template FOVs
        if [ ${FinalfMRIResolution} = "2" ] ; then
            ResampRefIm=$FSLDIR/data/standard/MNI152_T1_2mm
        elif [ ${FinalfMRIResolution} = "1" ] ; then
            ResampRefIm=$FSLDIR/data/standard/MNI152_T1_1mm
        else
            ${FSLDIR}/bin/flirt -interp spline -in ${T1wImage2} -ref ${T1wImage2} -applyisoxfm $FinalfMRIResolution -out ${onestep}/${T1wImageFile2}.${FinalfMRIResolution}
            ResampRefIm=${onestep}/${T1wImageFile2}.${FinalfMRIResolution} 
        fi
    
        if [ ! -e ${onestep}/${T1wImageFile2}.${FinalfMRIResolution}.nii.gz ] ; then
            echo "${i} Resting Functional One Step Resampling : Creating ${FinalfMRIResolution}mm standard space file for T1w image" 
            ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${T1wImage2} -r ${ResampRefIm} --premat=$FSLDIR/etc/flirtsch/ident.mat -o ${onestep}/${T1wImageFile2}.${FinalfMRIResolution}
            ${FSLDIR}/bin/imcp ${onestep}/${T1wImageFile2}.${FinalfMRIResolution} ${fMRIFolder}/${T1wImageFile2}.${FinalfMRIResolution}
        elif [ -e ${onestep}/${T1wImageFile2}.${FinalfMRIResolution}.nii.gz ] && [ -e ${fMRIFolder}/${T1wImageFile2}.${FinalfMRIResolution}.nii.gz ]; then
            echo "${i} Resting Functional One Step Resampling : ${FinalfMRIResolution}mm Standard Space File for T1w Image Created"
        else
            echo "${i} Resting Functional One Step Resampling : ERROR Creating ${FinalfMRIResolution}mm standard space file for T1w image" 
        fi
    
        if [ ! -e ${onestep}/${FreeSurferBrainMaskFile}.${FinalfMRIResolution}.nii.gz ] ; then
            echo "${i} Resting Functional One Step Resampling : Creating ${FinalfMRIResolution}mm standard space file for brain mask" 
            ${FSLDIR}/bin/applywarp --rel --interp=nn -i ${FreeSurferBrainMask}.nii.gz -r ${onestep}/${T1wImageFile2}.${FinalfMRIResolution} --premat=$FSLDIR/etc/flirtsch/ident.mat -o ${onestep}/${FreeSurferBrainMaskFile}.${FinalfMRIResolution}.nii.gz
            ${FSLDIR}/bin/imcp ${onestep}/${FreeSurferBrainMaskFile}.${FinalfMRIResolution} ${fMRIFolder}/${FreeSurferBrainMaskFile}.${FinalfMRIResolution}
        elif [ -e ${onestep}/${FreeSurferBrainMaskFile}.${FinalfMRIResolution}.nii.gz ] && [ -e ${fMRIFolder}/${FreeSurferBrainMaskFile}.${FinalfMRIResolution}.nii.gz ] ; then
            echo "${i} Resting Functional One Step Resampling : ${FinalfMRIResolution}mm Standard Space File for Brain Mask Created"
        else
            echo "${i} Resting Functional One Step Resampling : ERROR Creating ${FinalfMRIResolution}mm standard space file for brain mask" 
        fi
    
        #### Create versions of the biasfield (changing resolution) #### NEED TO WORK ON IT (SKIPPED BC THIS REQUIRES POSTFREESURFER)
        #if [ ! -e ${onestep}/${BiasFieldFile2}.${FinalfMRIResolution}.nii.gz ] ; then
        #    ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${BiasField2} -r ${onestep}/${FreeSurferBrainMaskFile}.${FinalfMRIResolution}.nii.gz --premat=$FSLDIR/etc/flirtsch/ident.mat -o ${onestep}/${BiasFieldFile2}.${FinalfMRIResolution}
        #    ${FSLDIR}/bin/fslmaths ${onestep}/${BiasFieldFile2}.${FinalfMRIResolution} -thr 0.1 ${onestep}/${BiasFieldFile2}.${FinalfMRIResolution}
        #    ${FSLDIR}/bin/imcp ${onestep}/${BiasFieldFile2}.${FinalfMRIResolution} ${fMRIFolder}/${BiasFieldFile2}.${FinalfMRIResolution}
        #fi
    
        #### Downsample warpfield (fMRI to standard) to increase speed ####
        if [ ! -e ${OutputTransform2}.nii.gz ] ; then 
            echo "${i} Resting Functional One Step Resampling : Downsampling REST-to-Standard warpfield" 
            ${FSLDIR}/bin/convertwarp --relout --rel --warp1=${OutputTransform} --warp2=${StructuralToStandard} --ref=${onestep}/${T1wImageFile2}.${FinalfMRIResolution} --out=${OutputTransform2}
        else
            echo "${i} Resting Functional One Step Resampling : REST-to-Standard Warpfield Downsampled" 
        fi
    
        ###Add stuff for RMS###
        if [ ! -e ${OutputInvTransform}.nii.gz ] ; then
            ${FSLDIR}/bin/invwarp -w ${OutputTransform2} -o ${OutputInvTransform} -r ${ScoutInputName}
        fi
        if [ ! -e ${ScoutInputName}_mask.nii.gz ] ; then
            ${FSLDIR}/bin/applywarp --rel --interp=nn -i ${FreeSurferBrainMask}.nii.gz -r ${ScoutInputName} -w ${OutputInvTransform} -o ${ScoutInputName}_mask.nii.gz
        fi

        #### Apply combined transformations to fMRI (combines gradient non-linearity distortion, motion correction, and registration to T1w space, but keeping fMRI resolution) ####
        ## using only mcflirt outputs for the time being (NEED TO EDIT FURHTER WITH ADDING MOTION CORRECTION TYPE AS OPTION INSTEAD OF RUNNING BOTH)
        if [ ! -d ${onestep}/prevols_mcflirt ] ; then
            mkdir -p ${onestep}/prevols_mcflirt
        fi
        if [ ! -d ${onestep}/postvols_mcflirt ] ; then
            mkdir -p ${onestep}/postvols_mcflirt
        fi

        if [ ! -e ${onestep}/postvols_mcflirt/vol0249_mask.nii.gz ]; then
            echo "${i} Resting Functional One Step Resampling : Applying combined transformations to fMRI (combines gradient non-linearity distortion, motion correction)"
            ${FSLDIR}/bin/fslsplit ${InputfMRI2} ${onestep}/prevols_mcflirt/vol -t
            FrameMergeSTRING=""
            FrameMergeSTRINGII=""
            for ((k=0; k < $NumFrames; k++)); do
                vnum=`${FSLDIR}/bin/zeropad $k 4`
                ###Add stuff for RMS###
                rmsdiff ${MotionMatrixFolder_mcflirt}/${MotionMatrixPrefix}${vnum} ${MotionMatrixFolder_mcflirt}/${MotionMatrixPrefix}0000 ${ScoutInputName} ${ScoutInputName}_mask.nii.gz | tail -n 1 >> ${fMRIFolder}/Movement_AbsoluteRMS.txt
                if [ $k -eq 0 ] ; then
                    echo "0" >> ${fMRIFolder}/Movement_RelativeRMS.txt
                else
                    rmsdiff ${MotionMatrixFolder_mcflirt}/${MotionMatrixPrefix}${vnum} $prevmatrix ${ScoutInputName} ${ScoutInputName}_mask.nii.gz | tail -n 1 >> ${fMRIFolder}/Movement_RelativeRMS.txt
                fi
                prevmatrix="${MotionMatrixFolder_mcflirt}/${MotionMatrixPrefix}${vnum}"
                ###Add stuff for RMS###
                # original script   ${FSLDIR}/bin/convertwarp --relout --rel --ref=${onestep}/prevols_mcflirt/vol${vnum}.nii.gz --warp1=${GradientDistortionField} --postmat=${MotionMatrixFolder_mcflirt}/${MotionMatrixPrefix}${vnum} --out=${MotionMatrixFolder_mcflirt}/${MotionMatrixPrefix}${vnum}_gdc_warp.nii.gz
                # original script   ${FSLDIR}/bin/convertwarp --relout --rel --ref=${onestep}/${T1wImageFile2}.${FinalfMRIResolution} --warp1=${MotionMatrixFolder_mcflirt}/${MotionMatrixPrefix}${vnum}_gdc_warp.nii.gz --warp2=${OutputTransform2} --out=${MotionMatrixFolder_mcflirt}/${MotionMatrixPrefix}${vnum}_all_warp.nii.gz
                # trial_1   convert_xfm -omat ${MotionMatrixFolder_mcflirt}/${MotionMatrixPrefix}${vnum}_motion2str.mat -concat ${DCFolder}/fMRI2str.mat ${MotionMatrixFolder_mcflirt}/${MotionMatrixPrefix}${vnum} 
                # trial_1   convert_xfm -omat ${MotionMatrixFolder_mcflirt}/${MotionMatrixPrefix}${vnum}_motion2std.mat -concat ${AtlasSpaceFolder}/xfms/acpc2MNILinear.mat ${MotionMatrixFolder_mcflirt}/${MotionMatrixPrefix}${vnum}_motion2str.mat 
                # trial_1   convertwarp --relout --rel --ref=${onestep}/${T1wImageFile2}.${FinalfMRIResolution} --premat=${MotionMatrixFolder_mcflirt}/${MotionMatrixPrefix}${vnum}_motion2std.mat --warp1=${OutputTransform2} --out=${MotionMatrixFolder_mcflirt}/${MotionMatrixPrefix}${vnum}_all_warp.nii.gz

                ${FSLDIR}/bin/convertwarp --relout --rel --ref=${onestep}/${T1wImageFile2}.${FinalfMRIResolution} --premat=${MotionMatrixFolder_mcflirt}/${MotionMatrixPrefix}${vnum} --warp1=${OutputTransform2} --out=${MotionMatrixFolder_mcflirt}/${MotionMatrixPrefix}${vnum}_all_warp.nii.gz


                ${FSLDIR}/bin/fslmaths ${onestep}/prevols_mcflirt/vol${vnum}.nii.gz -mul 0 -add 1 ${onestep}/prevols_mcflirt/vol${vnum}_mask.nii.gz     #(Q)BLANK IMAGE GENERATED
                ${FSLDIR}/bin/applywarp --rel --interp=spline --in=${onestep}/prevols_mcflirt/vol${vnum}.nii.gz --warp=${MotionMatrixFolder_mcflirt}/${MotionMatrixPrefix}${vnum}_all_warp.nii.gz --ref=${onestep}/${T1wImageFile2}.${FinalfMRIResolution} --out=${onestep}/postvols_mcflirt/vol${vnum}.nii.gz
                ${FSLDIR}/bin/applywarp --rel --interp=nn --in=${onestep}/prevols_mcflirt/vol${vnum}_mask.nii.gz --warp=${MotionMatrixFolder_mcflirt}/${MotionMatrixPrefix}${vnum}_all_warp.nii.gz --ref=${onestep}/${T1wImageFile2}.${FinalfMRIResolution} --out=${onestep}/postvols_mcflirt/vol${vnum}_mask.nii.gz
                FrameMergeSTRING="${FrameMergeSTRING}${onestep}/postvols_mcflirt/vol${vnum}.nii.gz " 
                FrameMergeSTRINGII="${FrameMergeSTRINGII}${onestep}/postvols_mcflirt/vol${vnum}_mask.nii.gz "
            done
        elif [  -e ${onestep}/postvols_mcflirt/vol0249_mask.nii.gz ]; then
            correct_count=500
            counting=`ls ${onestep}/postvols_mcflirt | wc -l`
            if [ $correct_count -eq $counting ] ; then
                echo "${i} Resting Functional One Step Resampling : Correctly applying combined transformations to fMRI (combines gradient non-linearity distortion, motion correction)"
            else
                echo "${i} Resting Functional One Step Resampling : ERROR applying combined transformations to fMRI (combines gradient non-linearity distortion, motion correction)"
            fi
        fi

        ##### Merge together results and restore the TR (saved beforehand) ####
        if [ ! -e ${OutputfMRI}.nii.gz ] ; then
            echo "${i} Resting Functional One Step Resampling : Merging together results and restoring the TR"
            ${FSLDIR}/bin/fslmerge -tr ${OutputfMRI} $FrameMergeSTRING $TR_vol
            ${FSLDIR}/bin/fslmerge -tr ${OutputfMRI}_mask $FrameMergeSTRINGII $TR_vol
            ${FSLDIR}/bin/fslmaths ${OutputfMRI}_mask -Tmin ${OutputfMRI}_mask
        elif [  -e ${OutputfMRI}.nii.gz ] && [ -e ${OutputfMRI}_mask.nii.gz ] ; then
            echo "${i} Resting Functional One Step Resampling : Results Merged and TR Restored"
        else
            echo "${i} Resting Functional One Step Resampling : ERROR Merging together results and restoring the TR"
        fi

        ##### Combine transformations: gradient non-linearity distortion + fMRI_dc to standard ####
        if [ ! -e ${onestep}/Scout_gdc_MNI_warp.nii.gz ] ; then
            echo "${i} Resting Fuctional One Step Resampling : Combining transformation (gradient non-linearity distortion + fMRI_dc to standard)"
            ${FSLDIR}/bin/convertwarp --relout --rel --ref=${onestep}/${T1wImageFile2}.${FinalfMRIResolution} --warp1=${GradientDistortionField} --warp2=${OutputTransform2} --out=${onestep}/Scout_gdc_MNI_warp.nii.gz
            ${FSLDIR}/bin/applywarp --rel --interp=spline --in=${ScoutInput} -w ${onestep}/Scout_gdc_MNI_warp.nii.gz -r ${onestep}/${T1wImageFile2}.${FinalfMRIResolution} -o ${ScoutOutput}
        elif [  -e ${onestep}/Scout_gdc_MNI_warp.nii.gz ] && [ -e ${ScoutOutput}.nii.gz ] ; then
            echo "${i} Resting Fuctional One Step Resampling : Transformation Combined (gradient non-linearity distortion + fMRI_dc to standard)"
        fi
    else
        #if 
        echo "${i} RESTING FUNCTIONAL ONE STEP RESAMPLING COMPLETED"
    fi

    ##### Create spline interpolated version of Jacobian  (T1w space, fMRI resolution) ####
    #OutputTransform is from gdc space to T1w space, ie, only fieldmap-based distortions (like topup)
    #output jacobian is both gdc and topup/fieldmap jacobian, but not the to MNI jacobian
    #JacobianIn was removed from inputs, now we just compute it from the combined warpfield of gdc and dc (NOT MNI)
    #compute combined warpfield, but don't use jacobian output because it has 8 frames for no apparent reason
    if [ ! -e ${onestep}/gdc_dc_jacobian.nii.gz ] ; then
        ${FSLDIR}/bin/convertwarp --relout --rel --ref=${OutputTransform} --warp1=${GradientDistortionField} --warp2=${OutputTransform} -o ${onestep}/gdc_dc_warp --jacobian=${onestep}/gdc_dc_jacobian
        ## but, convertwarp's jacobian is 8 frames - each combination of one-sided differences, so average them
        ${FSLDIR}/bin/fslmaths ${onestep}/gdc_dc_jacobian -Tmean ${onestep}/gdc_dc_jacobian
    fi
    if [ ! -e ${JacobianOut2}.nii.gz ] ; then
        ## resample it to MNI space
        ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${onestep}/gdc_dc_jacobian -r ${onestep}/${T1wImageFile2}.${FinalfMRIResolution} -w ${StructuralToStandard} -o ${JacobianOut2}
    fi

    ###Add stuff for RMS###
    if [ ! -e ${fMRIFolder}/Movement_RelativeRMS_mean.txt ] ; then
        cat ${fMRIFolder}/Movement_RelativeRMS.txt | awk '{ sum += $1} END { print sum / NR }' >> ${fMRIFolder}/Movement_RelativeRMS_mean.txt
    fi
    if [ ! -e ${fMRIFolder}/Movement_AbsoluteRMS_mean.txt ] ; then
        cat ${fMRIFolder}/Movement_AbsoluteRMS.txt | awk '{ sum += $1} END { print sum / NR }' >> ${fMRIFolder}/Movement_AbsoluteRMS_mean.txt
    fi

done









    
        



















