## Set diresctorpy paths
FSLDIR=/usr/local/fsl
HCPPIPEDIR=/Users/yoobin_kwak/preproc_BCS_001_002/Pipelines-master
CARET7DIR=/Users/yoobin_kwak/preproc_BCS_001_002/workbench/bin_macosx64

source $HCPPIPEDIR/global/scripts/log.shlib     # Logging related functions
source $HCPPIPEDIR/global/scripts/opts.shlib    # Command line option functions

for i in $@
do
    preproc_struct=${i}/Preprocess_Structure
    if [ ! -d ${preproc_struct} ]; then
        mkdir ${preproc_struct}
    fi
    pre_FS=${preproc_struct}/PreFreesurfer
    if [ ! -d ${pre_FS} ] ; then
        mkdir ${pre_FS}
    fi
    T1wFolder=${pre_FS}/T1w
    if [ ! -d ${T1wFolder} ] ; then
        mkdir ${T1wFolder}
        mkdir ${T1wFolder}/xfms
    fi
    T2wFolder=${pre_FS}/T2w
    if [ ! -d ${T2wFolder} ] ; then
        mkdir ${T2wFolder}
        mkdir ${T2wFolder}/xfms
    fi
    AtlasSpaceFolder=${pre_FS}/MNINonLinear
    if [ ! -d ${AtlasSpaceFolder} ] ; then
        mkdir ${AtlasSpaceFolder}
        mkdir ${AtlasSpaceFolder}/xfms
    fi


    T1wInputImage=${i}/T1*/20*
    T1wImage="T1w"
    T2wInputImage=${i}/T2*/20*
    T2wImage="T2w"
    T1wTemplate=$HCPPIPEDIR/global/templates/MNI152_T1_0.8mm.nii.gz
    T1wTemplateBrain=$HCPPIPEDIR/global/templates/MNI152_T1_0.8mm_brain.nii.gz
    T1wTemplate2mm=$HCPPIPEDIR/global/templates/MNI152_T1_2mm.nii.gz
    T2wTemplate=$HCPPIPEDIR/global/templates/MNI152_T2_0.8mm.nii.gz
    T2wTemplateBrain=$HCPPIPEDIR/global/templates/MNI152_T2_0.8mm_brain.nii.gz
    T2wTemplate2mm=$HCPPIPEDIR/global/templates/MNI152_T2_2mm.nii.gz
    TemplateMask=${HCPPIPEDIR}/global/templates/MNI152_T1_0.8mm_brain_mask.nii.gz
    Template2mmMask=${HCPPIPEDIR}/global/templates/MNI152_T1_2mm_brain_mask_dil.nii.gz
    FNIRTConfig=$FSLDIR/etc/flirtsch/T1_2_MNI152_2mm.cnf
    
    Modalities="T1w T2w"
    for TXw in ${Modalities}
    do    
        # set up appropriate input variables
        if [ $TXw = T1w ] ; then
            TXwInputImage="${T1wInputImage}"
            TXwFolder=${T1wFolder}
            TXwImage=${T1wImage}
            TXwTemplate=${T1wTemplate}
            TXwTemplate2mm=${T1wTemplate2mm}
        else
            TXwInputImage="${T2wInputImage}"
            TXwFolder=${T2wFolder}
            TXwImage=${T2wImage}
            TXwTemplate=${T2wTemplate}
            TXwTemplate2mm=${T2wTemplate2mm}
        fi
        OutputTXwImageSTRING=""
        
        GradientDistortionCoeffs=NONE       ## File containing gradient distortion coefficients, Set to "NONE" to turn off (we currently don't have the info)
        echo "${i} PRE_FREESURFER : NOT PERFORMING GRADIENT DISTORTION CORRECTION"
       
        if [ ! -f ${TXwFolder}/${TXwImage}1_gdc.nii.gz ] ; then
            j=1
            for Image in $TXwInputImage
            do
                ${FSLDIR}/bin/fslreorient2std $Image ${TXwFolder}/${TXwImage}${j}_gdc   
                OutputTXwImageSTRING="${OutputTXwImageSTRING}${TXwFolder}/${TXwImage}${j}_gdc "
                j=$(($j+1))
            done
        fi

        echo "${i} PRE_FREESURFER: Averaging neither T1 nor T2"
        if [ ! -f ${TXwFolder}/${TXwImage}.nii.gz ] ; then
            ${FSLDIR}/bin/imcp ${TXwFolder}/${TXwImage}1_gdc ${TXwFolder}/${TXwImage}
        fi

#### ACPC align T1w or T2w image to 0.7mm MNI Template to create native volume space ####
        wd_acpc=${TXwFolder}/ACPCAlignment
        if [ ! -d ${wd_acpc} ] ; then
            mkdir ${wd_acpc}
        fi
        if [ ! -f ${TXwFolder}/${TXwImage}_acpc.nii.gz ] ; then
            echo "${i} PRE_FREESURFER : Aligning ${TXw} image to 0.7mm MNI ${TXw}Template to create native volume space"
            # Crop the FOV
            ${FSLDIR}/bin/robustfov -i ${TXwFolder}/${TXwImage} -m ${wd_acpc}/roi2full.mat -r ${wd_acpc}/robustroi.nii.gz
            # Invert the matrix (to get full FOV to ROI)
            ${FSLDIR}/bin/convert_xfm -omat ${wd_acpc}/full2roi.mat -inverse ${wd_acpc}/roi2full.mat
            # Register cropped image to MNI152 (12 DOF)
            ${FSLDIR}/bin/flirt -interp spline -in ${wd_acpc}/robustroi.nii.gz -ref ${TXwTemplate} -omat ${wd_acpc}/roi2std.mat -out ${wd_acpc}/acpc_final.nii.gz -searchrx -30 30 -searchry -30 30 -searchrz -30 30
            # Concatenate matrices to get full FOV to MNI
            ${FSLDIR}/bin/convert_xfm -omat ${wd_acpc}/full2std.mat -concat ${wd_acpc}/roi2std.mat ${wd_acpc}/full2roi.mat
            # Get a 6 DOF approximation which does the ACPC alignment (AC, ACPC line, and hemispheric plane)
            ${FSLDIR}/bin/aff2rigid ${wd_acpc}/full2std.mat ${TXwFolder}/xfms/acpc.mat
            # Create a resampled image (ACPC aligned) using spline interpolation
            ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${TXwFolder}/${TXwImage} -r ${TXwTemplate} --premat=${TXwFolder}/xfms/acpc.mat -o ${TXwFolder}/${TXwImage}_acpc 
        else
            echo "${i} PRE_FREESURFER : ${TXw} ACPC aligned"
        fi

#### Brain Extraction(FNIRT-based Masking) ####
        wd_bet=${TXwFolder}/BrainExtraction_FNIRTbased
        if [ ! -d ${wd_bet} ] ; then
            mkdir ${wd_bet}
        fi
        if [ ! -f ${TXwFolder}/${TXwImage}_acpc_brain.nii.gz ] ; then
            echo "${i} PRE_FREESURFER : Performing Brain Extraction using FNIRT-based Masking"
            BaseName=`${FSLDIR}/bin/remove_ext ${TXwFolder}/${TXwImage}_acpc`;
            BaseName=`basename $BaseName`;
            # Register to 2mm reference image (linear then non-linear)
            ${FSLDIR}/bin/flirt -interp spline -dof 12 -in ${TXwFolder}/${TXwImage}_acpc -ref ${TXwTemplate2mm} -omat ${wd_bet}/roughlin.mat -out ${wd_bet}/"$BaseName"_to_MNI_roughlin.nii.gz -nosearch
            ${FSLDIR}/bin/fnirt --in=${TXwFolder}/${TXwImage}_acpc --ref=${T1wTemplate2mm} --aff=${wd_bet}/roughlin.mat --refmask=${Template2mmMask} --fout=${wd_bet}/str2standard.nii.gz --jout=${wd_bet}/NonlinearRegJacobians.nii.gz --refout=${wd_bet}/IntensityModulatedT1.nii.gz --iout=${wd_bet}/"$BaseName"_to_MNI_nonlin.nii.gz --logout=${wd_bet}/NonlinearReg.txt --intout=${wd_bet}/NonlinearIntensities.nii.gz --cout=${wd_bet}/NonlinearReg.nii.gz --config=$FNIRTConfig
                # Overwrite the image output from FNIRT with a spline interpolated highres version
            ${FSLDIR}/bin/applywarp --rel --interp=spline --in=${TXwFolder}/${TXwImage}_acpc --ref=${TXwTemplate} -w ${wd_bet}/str2standard.nii.gz --out=${wd_bet}/"$BaseName"_to_MNI_nonlin.nii.gz
            # Invert warp and transform dilated brain mask back into native space, and use it to mask input image
            # {TXwFolder}/${TXwImage}_acpc and reference spaces are the same, using 2mm reference to save time
            ${FSLDIR}/bin/invwarp --ref=${TXwTemplate2mm} -w ${wd_bet}/str2standard.nii.gz -o ${wd_bet}/standard2str.nii.gz
            ${FSLDIR}/bin/applywarp --rel --interp=nn --in=${TemplateMask} --ref=${TXwFolder}/${TXwImage}_acpc -w ${wd_bet}/standard2str.nii.gz -o ${TXwFolder}/${TXwImage}_acpc_brain_mask
            ${FSLDIR}/bin/fslmaths ${TXwFolder}/${TXwImage}_acpc -mas ${TXwFolder}/${TXwImage}_acpc_brain_mask ${TXwFolder}/${TXwImage}_acpc_brain 
        else
            echo "${i} PRE_FREESURFER : ${TXw} Brain Extraced using FNIRT"
        fi
    done
        
#### T2w to T1w Registration (Optional Readout Distortion Correction not included in the current script) ####
        wd_T2T1Reg=${T2wFolder}/T2wToT1wReg        
        if [ ! -d ${wd_T2T1Reg} ] ;then
            mkdir ${wd_T2T1Reg}
        fi
        if [ ! -f ${T1wFolder}/${T2wImage}_acpc_dc.nii.gz ] ; then
            echo "${i} PRE_FREESURFER : Performing T2w to T1w Registration"      ## Not performing readout distortion correction 
            t1=${T1wFolder}/${T1wImage}_acpc 
            t1brain=${T1wFolder}/${T1wImage}_acpc_brain 
            t2=${T2wFolder}/${T2wImage}_acpc 
            t2brain=${T2wFolder}/${T2wImage}_acpc_brain 
            ot1=${T1wFolder}/${T1wImage}_acpc_dc 
            ot1brain=${T1wFolder}/${T1wImage}_acpc_dc_brain 
            ot1warp=${T1wFolder}/xfms/${T1wImage}_dc 
            ot2=${T1wFolder}/${T2wImage}_acpc_dc 
            ot2warp=${T1wFolder}/xfms/${T2wImage}_reg_dc
            t1brainFile=`basename "$t1brain"`

            cp "$t1brain".nii.gz ${wd_T2T1Reg}/"$t1brainFile".nii.gz
            ${FSLDIR}/bin/epi_reg --epi="$t2brain" --t1="$t1" --t1brain=${wd_T2T1Reg}/"$t1brainFile" --out=${wd_T2T1Reg}/T2w2T1w
            ${FSLDIR}/bin/applywarp --rel --interp=spline --in="$t2" --ref="$t1" --premat=${wd_T2T1Reg}/T2w2T1w.mat --out=${wd_T2T1Reg}/T2w2T1w
            ${FSLDIR}/bin/fslmaths ${wd_T2T1Reg}/T2w2T1w -add 1 ${wd_T2T1Reg}/T2w2T1w -odt float
            cp "$t1".nii.gz "$ot1".nii.gz
            cp "$t1brain".nii.gz "$ot1brain".nii.gz
            ${FSLDIR}/bin/fslmerge -t $ot1warp "$t1".nii.gz "$t1".nii.gz "$t1".nii.gz
            ${FSLDIR}/bin/fslmaths $ot1warp -mul 0 $ot1warp
            cp ${wd_T2T1Reg}/T2w2T1w.nii.gz "$ot2".nii.gz
         ${FSLDIR}/bin/convertwarp --relout --rel -r "$ot2".nii.gz -w $ot1warp --postmat=${wd_T2T1Reg}/T2w2T1w.mat --out="$ot2warp"
        else
            echo "${i} PRE_FREESURFER : T2w Registered to T1w (without readout distortion correction but output names as if)"
        fi

#### Bias Field Correction: Calculate bias field using square root of the product of T1w and T2w iamges ####
        wd_bfc=${T1wFolder}/BiasFieldCorrection_sqrtT1wXT1w 
        if [ ! -d ${wd_bfc} ] ; then
            mkdir ${wd_bfc}
        fi
        if [ ! -f ${T1wFolder}/${T2wImage}_acpc_dc_restore_brain.nii.gz ] ;then
            echo "${i} PRE_FREESURFER : Performing Bias Field Correction"
            T1im=${T1wFolder}/${T1wImage}_acpc_dc 
            T1brain=${T1wFolder}/${T1wImage}_acpc_dc_brain 
            T2im=${T1wFolder}/${T2wImage}_acpc_dc
            obias=${T1wFolder}/BiasField_acpc_dc            ## output bias field image 
            oT1im=${T1wFolder}/${T1wImage}_acpc_dc_restore          ## output corrected T1 image
            oT1brain=${T1wFolder}/${T1wImage}_acpc_dc_restore_brain             ## output corrected T1 brain image
            oT2im=${T1wFolder}/${T2wImage}_acpc_dc_restore          ## output corrected T2 image      
            oT2brain=${T1wFolder}/${T2wImage}_acpc_dc_restore_brain             ## output corrected T2 brain
            BiasFieldSmoothingSigma=5       ## default
            Factor=0.5      ## default
        
            # Form sqrt(T1w*T2w), mask this and normalise by the mean                                                                 
            ${FSLDIR}/bin/fslmaths $T1im -mul $T2im -abs -sqrt ${wd_bfc}/T1wmulT2w.nii.gz -odt float
            ${FSLDIR}/bin/fslmaths ${wd_bfc}/T1wmulT2w.nii.gz -mas $T1brain ${wd_bfc}/T1wmulT2w_brain.nii.gz
            meanbrainval=`${FSLDIR}/bin/fslstats ${wd_bfc}/T1wmulT2w_brain.nii.gz -M`
            ${FSLDIR}/bin/fslmaths ${wd_bfc}/T1wmulT2w_brain.nii.gz -div $meanbrainval ${wd_bfc}/T1wmulT2w_brain_norm.nii.gz
    
            # Smooth the normalised sqrt image, using within-mask smoothing : s(Mask*X)/s(Mask)
            ${FSLDIR}/bin/fslmaths ${wd_bfc}/T1wmulT2w_brain_norm.nii.gz -bin -s $BiasFieldSmoothingSigma ${wd_bfc}/SmoothNorm_s${BiasFieldSmoothingSigma}.nii.gz
            ${FSLDIR}/bin/fslmaths ${wd_bfc}/T1wmulT2w_brain_norm.nii.gz -s $BiasFieldSmoothingSigma -div ${wd_bfc}/SmoothNorm_s${BiasFieldSmoothingSigma}.nii.gz ${wd_bfc}/T1wmulT2w_brain_norm_s${BiasFieldSmoothingSigma}.nii.gz
          
            # Divide normalised sqrt image by smoothed version (to do simple bias correction)
            ${FSLDIR}/bin/fslmaths ${wd_bfc}/T1wmulT2w_brain_norm.nii.gz -div ${wd_bfc}/T1wmulT2w_brain_norm_s$BiasFieldSmoothingSigma.nii.gz ${wd_bfc}/T1wmulT2w_brain_norm_modulate.nii.gz
      
            # Create a mask using a threshold at Mean - 0.5*Stddev, with filling of holes to remove any non-grey/white tissue.
            STD=`${FSLDIR}/bin/fslstats ${wd_bfc}/T1wmulT2w_brain_norm_modulate.nii.gz -S`
            echo $STD
            MEAN=`${FSLDIR}/bin/fslstats ${wd_bfc}/T1wmulT2w_brain_norm_modulate.nii.gz -M`
            echo $MEAN
            Lower=`echo "$MEAN - ($STD * $Factor)" | bc -l`
            echo $Lower
            ${FSLDIR}/bin/fslmaths ${wd_bfc}/T1wmulT2w_brain_norm_modulate -thr $Lower -bin -ero -mul 255 ${wd_bfc}/T1wmulT2w_brain_norm_modulate_mask
            ${CARET7DIR}/wb_command -volume-remove-islands ${wd_bfc}/T1wmulT2w_brain_norm_modulate_mask.nii.gz ${wd_bfc}/T1wmulT2w_brain_norm_modulate_mask.nii.gz
                  
            # Extrapolate normalised sqrt image from mask region out to whole FOV
            ${FSLDIR}/bin/fslmaths ${wd_bfc}/T1wmulT2w_brain_norm.nii.gz -mas ${wd_bfc}/T1wmulT2w_brain_norm_modulate_mask.nii.gz -dilall ${wd_bfc}/bias_raw.nii.gz -odt float
            ${FSLDIR}/bin/fslmaths ${wd_bfc}/bias_raw.nii.gz -s $BiasFieldSmoothingSigma $obias
                
            # Use bias field output to create corrected images
            ${FSLDIR}/bin/fslmaths $T1im -div $obias -mas $T1brain $oT1brain -odt float
            ${FSLDIR}/bin/fslmaths $T1im -div $obias $oT1im -odt float
            ${FSLDIR}/bin/fslmaths $T2im -div $obias -mas $T1brain $oT2brain -odt float
            ${FSLDIR}/bin/fslmaths $T2im -div $obias $oT2im -odt float
        else
            echo "${i} PRE_FREESURFER : Bias Field Corrected "
        fi

#### Atlas Registration to MNI152: FLIRT + FNIRT Also applies registration to T1w and T2w images ####
        if [ ! -f ${AtlasSpaceFolder}/${T2wImage}_restore_brain.nii.gz ] ; then
            echo "${i} PRE_FREESURFER : Performing Atlas Registration to MNI152 (FLIRT and FNIRT)"
            t1=${T1wFolder}/${T1wImage}_acpc_dc 
            t1rest=${T1wFolder}/${T1wImage}_acpc_dc_restore 
            t1restbrain=${T1wFolder}/${T1wImage}_acpc_dc_restore_brain 
            t2=${T1wFolder}/${T2wImage}_acpc_dc 
            t2rest=${T1wFolder}/${T2wImage}_acpc_dc_restore 
            t2restbrain=${T1wFolder}/${T2wImage}_acpc_dc_restore_brain 
            owarp=${AtlasSpaceFolder}/xfms/acpc_dc2standard.nii.gz 
            oinvwarp=${AtlasSpaceFolder}/xfms/standard2acpc_dc.nii.gz 
            ot1=${AtlasSpaceFolder}/${T1wImage} 
            ot1rest=${AtlasSpaceFolder}/${T1wImage}_restore 
            ot1restbrain=${AtlasSpaceFolder}/${T1wImage}_restore_brain 
            ot2=${AtlasSpaceFolder}/${T2wImage}
            ot2rest=${AtlasSpaceFolder}/${T2wImage}_restore 
            ot2restbrain=${AtlasSpaceFolder}/${T2wImage}_restore_brain 

            t1restBasename=`remove_ext $t1rest`;
            t1restBasename=`basename $t1restBasename`;
            t1restbrainBasename=`remove_ext $t1restbrain`;
            t1restbrainBasename=`basename $t1restbrainBasename`;
        
            # Linear then non-linear registration to MNI
            ${FSLDIR}/bin/flirt -interp spline -dof 12 -in ${t1restbrain} -ref ${T1wTemplateBrain} -omat ${AtlasSpaceFolder}/xfms/acpc2MNILinear.mat -out ${AtlasSpaceFolder}/xfms/${t1restbrainBasename}_to_MNILinear
            ${FSLDIR}/bin/fnirt --in=${t1rest} --ref=${T1wTemplate2mm} --aff=${AtlasSpaceFolder}/xfms/acpc2MNILinear.mat --refmask=${Template2mmMask} --fout=${owarp} --jout=${AtlasSpaceFolder}/xfms/NonlinearRegJacobians.nii.gz --refout=${AtlasSpaceFolder}/xfms/IntensityModulatedT1.nii.gz --iout=${AtlasSpaceFolder}/xfms/2mmReg.nii.gz --logout=${AtlasSpaceFolder}/xfms/NonlinearReg.txt --intout=${AtlasSpaceFolder}/xfms/NonlinearIntensities.nii.gz --cout=${AtlasSpaceFolder}/xfms/NonlinearReg.nii.gz --config=${FNIRTConfig}
  
            # Input and reference spaces are the same, using 2mm reference to save time
            ${FSLDIR}/bin/invwarp -w ${owarp} -o ${oinvwarp} -r ${T1wTemplate2mm}
  
            # T1w set of warped outputs (brain/whole-head + restored/orig)
            ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${t1} -r ${T1wTemplate} -w ${owarp} -o ${ot1}
            ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${t1rest} -r ${T1wTemplate} -w ${owarp} -o ${ot1rest}
            ${FSLDIR}/bin/applywarp --rel --interp=nn -i ${t1restbrain} -r ${T1wTemplate} -w ${owarp} -o ${ot1restbrain}
            ${FSLDIR}/bin/fslmaths ${ot1rest} -mas ${ot1restbrain} ${ot1restbrain}
 
            # T2w set of warped outputs (brain/whole-head + restored/orig)
            ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${t2} -r ${T1wTemplate} -w ${owarp} -o ${ot2}
            ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${t2rest} -r ${T1wTemplate} -w ${owarp} -o ${ot2rest}
            ${FSLDIR}/bin/applywarp --rel --interp=nn -i ${t2restbrain} -r ${T1wTemplate} -w ${owarp} -o ${ot2restbrain}
            ${FSLDIR}/bin/fslmaths ${ot2rest} -mas ${ot2restbrain} ${ot2restbrain}
        else
            echo "${i} PRE_FREESURFER : T1w and T2w Registered to MNI"
        fi

        echo "${i} PREEFREESURFER PROCESSING COMPLETED"

done


