HCPPIPEDIR=Pipelines-master
HCPPIPEDIR_Global=${HCPPIPEDIR}/global/scripts
FSLDIR=/usr/share/fsl/5.0

for i in $@
do
    preproc_diff=${i}/Preprocess_Diffusion
    if [ ! -d ${preproc_diff} ]
    then
        mkdir ${preproc_diff}
    fi
    rawdir=${preproc_diff}/raw
    if [ ! -d ${rawdir} ]
    then
        mkdir ${rawdir}
    fi
    if [ ! -f ${rawdir}/LR_B3000.bvec ]
    then
        cp -r ${i}/DTI_BLIP_RL*/*.nii.gz ${rawdir}/RL_BLIP.nii.gz 
        cp -r ${i}/DTI_BLIP_RL*/*.bval ${rawdir}/RL_BLIP.bval 
        cp -r ${i}/DTI_BLIP_RL*/*.bvec ${rawdir}/RL_BLIP.bvec 
        cp -r ${i}/DTI_BLIP_LR*/*.nii.gz ${rawdir}/LR_BLIP.nii.gz 
        cp -r ${i}/DTI_BLIP_LR*/*.bval ${rawdir}/LR_BLIP.bval 
        cp -r ${i}/DTI_BLIP_LR*/*.bvec ${rawdir}/LR_BLIP.bvec
        cp -r ${i}/DTI_MB3_LR_B1000*/*.nii.gz ${rawdir}/LR_B1000.nii.gz
        cp -r ${i}/DTI_MB3_LR_B1000*/*.bval ${rawdir}/LR_B1000.bval
        cp -r ${i}/DTI_MB3_LR_B1000*/*.bvec ${rawdir}/LR_B1000.bvec
        cp -r ${i}/DTI_MB3_LR_B2000*/*.nii.gz ${rawdir}/LR_B2000.nii.gz
        cp -r ${i}/DTI_MB3_LR_B2000*/*.bval ${rawdir}/LR_B2000.bval
        cp -r ${i}/DTI_MB3_LR_B2000*/*.bvec ${rawdir}/LR_B2000.bvec
        cp -r ${i}/DTI_MB3_LR_B3000*/*.nii.gz ${rawdir}/LR_B3000.nii.gz
        cp -r ${i}/DTI_MB3_LR_B3000*/*.bval ${rawdir}/LR_B3000.bval
        cp -r ${i}/DTI_MB3_LR_B3000*/*.bvec ${rawdir}/LR_B3000.bvec
    else
        echo ${i} raw diffusion directory exists
    fi

    topupdir=${preproc_diff}/topup
    if [ ! -d ${topupdir} ]
    then
        mkdir ${topupdir}
    fi
    eddydir=${preproc_diff}/eddy
    if [ ! -d ${eddydir} ]
    then
        mkdir ${eddydir}
    fi

    b0maxbval=50        ## from DiffPreprocPipeline.sh
    DEFAULT_DEGREES_OF_FREEDOM=6        ## from DiffPreprocPipeline.sh
    b0dist=45       ## from DiffPreprocPipeline_PreEddy.sh
    echo_spacing=0.87       ## in msec
    PEdir=-x        ## RL
        ## basePos=RL
        ## baseNeg=LR
    grappa_factor=3

#### Compute Total_readout in secs with up to 6 decimal places (basic_preproc.sh) ####
    dimP=`${FSLDIR}/bin/fslval ${rawdir}/RL_BLIP.nii.gz dim1`
    nPEsteps=$(($dimP - 1))
    ro_time=`echo "${echo_spacing} * (${nPEsteps}/${grappa_factor})" | bc -l`
    ro_time=`echo "scale=6; ${ro_time} / 1000" | bc -l`
    echo "basic_preproc: Total readout time is $ro_time secs"

#### Rescale series to ensure consistency across baseline intensities (basic_preproc.sh) ####
    entry_cnt=0
    for entry in ${rawdir}/RL_BLIP.nii.gz ${rawdir}/LR_BLIP.nii.gz ${rawdir}/LR_B1000.nii.gz ${rawdir}/LR_B2000.nii.gz ${rawdir}/LR_B3000.nii.gz 
    do
        basename=`imglob ${entry}`
        if [ ! -f ${basename}_mean.nii.gz ]
        then
            ${FSLDIR}/bin/fslmaths ${entry} -Xmean -Ymean -Zmean ${basename}_mean
        else
            echo ${i} ${basename} diffusion mean file created
        fi
        Posbvals=`cat ${basename}.bval`
#        echo "${basename} Posbvals: ${Posbvals}"
        
        #### extract all b0s for the series ####
        #if [ ! -f ${rawdir}/RL_BLIP_b0_0000.nii.gz ] && [ ! -f ${rawdir}/LR_BLIP_b0_0000.nii.gz ]
        if [ ! -f ${basename}_b0_0000.nii.gz ]
        then
            mcnt=0
            for p in ${Posbvals}
            do
                echo "Posbvals p: ${p}"
                cnt=`$FSLDIR/bin/zeropad $mcnt 4`
                echo "cnt: ${cnt}"
                if [ $p -lt ${b0maxbval} ]
                then
                    fslroi ${basename}_mean ${basename}_b0_${cnt} ${mcnt} 1
                fi
                mcnt=$((${mcnt} + 1))
            done
        else 
            echo ${i} ${basename} all b0 extracted for the series
        fi
        if [ ! -f ${basename}_merged_mean.nii.gz ]
        then
            ${FSLDIR}/bin/fslmerge -t ${basename}_merged `echo ${basename}_b0_????.nii*`        
            ${FSLDIR}/bin/fslmaths ${basename}_merged -Tmean ${basename}_merged_mean       ## This is the mean baseline b0 intensity for the series
            if [ ${entry_cnt} -eq 0 ]
            then
                rescale=`fslmeants -i ${basename}_merged_mean`
            else
                scaleS=`fslmeants -i ${basename}_merged_mean`
                ${FSLDIR}/bin/fslmaths ${basename} -mul ${rescale} -div ${scaleS} ${basename}_new
#               ${FSLDIR}/bin/imrm ${basename}
#               ${FSLDIR}/bin/immv ${basename}_new ${basename}
            fi
            entry_cnt=$((${entry_cnt} + 1))
        else
            echo ${i} ${basename} mean baseline b0 intensieity for the series processes
        fi
    done

#### Extract b0 and write text files for Topup and Eddy (basic_preproc.sh) ####
    scount=1
    indcount=0
    for entry in ${rawdir}/RL_BLIP.nii.gz
    do
        basename=`imglob ${entry}`
        Posbvals=`cat ${basename}.bval`
        count=0
        count3=$((${b0dist} + 1))
        for p in ${Posbvals}
        do
            if [ ! -f ${rawdir}/Pos_b0_0000.nii.gz ] && [ ! -f ${rawdir}/acqparams.txt ] 
            then
                if [ $p -lt ${b0maxbval} ] && [ ${count3} -gt ${b0dist} ]
                then
                    cnt=`$FSLDIR/bin/zeropad $indcount 4`
                    echo "Extracting Pos Volume $count from ${entry} as a b=0. Measured b=$p" >>${rawdir}/extractedb0.txt
                    $FSLDIR/bin/fslroi ${entry} ${rawdir}/Pos_b0_${cnt} ${count} 1
                    echo 1 0 0 ${ro_time} >> ${rawdir}/acqparams.txt
                    indcount=$((${indcount} + 1))
                    count3=0
                fi
#                if [ ! -f ${rawdir}/index.txt ]        ## using series_indext.txt instead of indext.txt for Eddy
#                then
#                    echo ${indcount} >> ${rawdir}/index.txt
#                    count3=$((${count3} + 1))
#                else
#                    echo ${i} ${basename} index.txt processed
#                fi
            else
                echo ${i} ${basename} ${p} b0 extracted and acqparams.txt processed
            fi
            count=$((${count} + 1))
        done

        if [ ! -f ${rawdir}/series_index.txt ]
        then
            sesdimt=`${FSLDIR}/bin/fslval ${entry} dim4`
            for (( j=0; j<${sesdimt}; j++ ))
            do
                echo ${scount} >> ${rawdir}/series_index.txt
            done
            scount=$((${scount} + 1))
        else
            echo ${i} ${basename} series_index.txt processsed
        fi
    done

    Poscount=${indcount}
    indcount=0
    for entry in ${rawdir}/LR_BLIP.nii.gz ${rawdir}/LR_B1000_new.nii.gz ${rawdir}/LR_B2000_new.nii.gz ${rawdir}/LR_B3000_new.nii.gz 
    do
        basename=`imglob ${entry}`
        edit_basename=`echo ${basename} | sed 's/_new//'`
        Negbvals=`cat ${edit_basename}.bval`
        count=0
        count3=$((${b0dist} + 1))
        for n in ${Negbvals}
        do
            correct_inputs=5
            check_acqparams=`cat ${rawdir}/acqparams.txt | wc -l` 
            check_extractedb0=`cat ${rawdir}/extractedb0.txt |wc -l`
            if [ $check_acqparams -eq $correct_inputs ] && [ $check_extractedb0 -eq $correct_inputs ]
            then
                echo ${i} b0 extraction and acqparam.txt completed : to be used for topup/eddy
            elif [ $check_acqparams -gt $correct_inputs ] && [ $check_extractedb0 -gt $correct_inputs ]
            then
                echo ERROR : ${i} need to re-proecess acqparams.txt and extractedb0.txt
            elif [ $check_acqparams -lt $correct_inputs ] && [ $check_extractedb0 -lt $correct_inputs ]
            then
                if [ $n -lt ${b0maxbval} ] && [ ${count3} -gt ${b0dist} ]
                then
                    cnt=`$FSLDIR/bin/zeropad $indcount 4`
                    echo "Extracting Neg Volume $count from ${entry} as a b=0. Measured b=$n" >>${rawdir}/extractedb0.txt
                    $FSLDIR/bin/fslroi ${entry} ${rawdir}/Neg_b0_${cnt} ${count} 1
                    echo -1 0 0 ${ro_time} >> ${rawdir}/acqparams.txt
                    indcount=$((${indcount} + 1))
                    count3=0
                fi
#                echo $((${indcount} + ${Poscount})) >>${rawdir}/index.txt
#                count3=$((${count3} + 1))
            fi
            count=$((${count} + 1))
        done

        correct_index=131
        check_index=`cat ${rawdir}/series_index.txt | wc -l`
        if [ $correct_index -eq $check_index ]
        then
            echo ${i} series_index.txt completed : to be used for eddy
        elif [ $correct_index -gt $check_index ]
        then
            echo ERROR : ${i} need to re-process sereis_index.txt 
        elif [ $correct_index -lt $check_index ]
        then
            sesdimt=`${FSLDIR}/bin/fslval ${entry} dim4`        ## this part needs to be "improved"
            for (( j=0; j<${sesdimt}; j++ ))
            do                
                echo ${scount} >> ${rawdir}/series_index.txt
            done
            scount=$((${scount} + 1))
        fi
    done

#### Merge files and correct number of slices (basic_preproc.sh) ####    
    if [ ! -f ${rawdir}/Pos_Neg.bvecs ] || [ ! -f ${eddydir}/Pos_Neg.bvecs ]
    then
        ${FSLDIR}/bin/fslmerge -t ${rawdir}/Pos_b0 `${FSLDIR}/bin/imglob ${rawdir}/Pos_b0_????.*`
        ${FSLDIR}/bin/fslmerge -t ${rawdir}/Neg_b0 `${FSLDIR}/bin/imglob ${rawdir}/Neg_b0_????.*`
        ${FSLDIR}/bin/fslmerge -t ${rawdir}/Pos ${rawdir}/RL_BLIP.nii.gz
        ${FSLDIR}/bin/fslmerge -t ${rawdir}/Neg ${rawdir}/LR_BLIP.nii.gz ${rawdir}/LR_B1000.nii.gz ${rawdir}/LR_B2000.nii.gz ${rawdir}/LR_B3000.nii.gz

        paste ${rawdir}/RL_BLIP.bval >${rawdir}/Pos.bval
        paste ${rawdir}/RL_BLIP.bvec >${rawdir}/Pos.bvec
        paste ${rawdir}/LR_BLIP.bval ${rawdir}/LR_B1000.bval ${rawdir}/LR_B2000.bval ${rawdir}/LR_B3000.bval >${rawdir}/Neg.bval
        paste ${rawdir}/LR_BLIP.bvec ${rawdir}/LR_B1000.bvec ${rawdir}/LR_B2000.bvec ${rawdir}/LR_B3000.bvec >${rawdir}/Neg.bvec

#       dimz=`${FSLDIR}/bin/fslval ${rawdir}/Pos dim3`         lines 209~224 error, but dimz=66 for BCS so skipped
#       if [ `isodd $dimz` -eq 1 ]
#       then
#           ${FSLDIR}/bin/fslroi ${rawdir}/Pos ${rawdir}/Posn 0 -1 0 -1 1 -1
#           ${FSLDIR}/bin/fslroi ${rawdir}/Neg ${rawdir}/Negn 0 -1 0 -1 1 -1
#           ${FSLDIR}/bin/fslroi ${rawdir}/Pos_b0 ${rawdir}/Pos_b0n 0 -1 0 -1 1 -1
#           ${FSLDIR}/bin/fslroi ${rawdir}/Neg_b0 ${rawdir}/Neg_b0n 0 -1 0 -1 1 -1
#           ${FSLDIR}/bin/imrm ${rawdir}/Pos
#           ${FSLDIR}/bin/imrm ${rawdir}/Neg
#           ${FSLDIR}/bin/imrm ${rawdir}/Pos_b0
#           ${FSLDIR}/bin/imrm ${rawdir}/Neg_b0
#           ${FSLDIR}/bin/immv ${rawdir}/Posn ${rawdir}/Pos
#           ${FSLDIR}/bin/immv ${rawdir}/Negn ${rawdir}/Neg
#           ${FSLDIR}/bin/immv ${rawdir}/Pos_b0n ${rawdir}/Pos_b0
#           ${FSLDIR}/bin/immv ${rawdir}/Neg_b0n ${rawdir}/Neg_b0
#       fi
        ${FSLDIR}/bin/fslmerge -t ${rawdir}/Pos_Neg_b0 ${rawdir}/Pos_b0 ${rawdir}/Neg_b0
        ${FSLDIR}/bin/fslmerge -t ${rawdir}/Pos_Neg ${rawdir}/Pos ${rawdir}/Neg
        paste ${rawdir}/Pos.bval ${rawdir}/Neg.bval >${rawdir}/Pos_Neg.bvals
        paste ${rawdir}/Pos.bvec ${rawdir}/Neg.bvec >${rawdir}/Pos_Neg.bvecs
    else
        echo ${i} diffusion data basic preprocessing completed
    fi

#### Move files to appropriate directories (basic_preproc.sh) ####
    if [ -f ${rawdir}/Neg.bvec ] &&  [ ! -f ${eddydir}/Neg.bvec ]
    then
        cp ${rawdir}/extractedb0.txt ${topupdir}/extractedb0.txt
        cp ${rawdir}/acqparams.txt ${topupdir}/acqparams.txt
        ${FSLDIR}/bin/immv ${rawdir}/Pos_Neg_b0 ${topupdir}
        ${FSLDIR}/bin/immv ${rawdir}/Pos_b0 ${topupdir}
        ${FSLDIR}/bin/immv ${rawdir}/Neg_b0 ${topupdir}

        cp ${topupdir}/acqparams.txt ${eddydir}
#        mv ${rawdir}/index.txt ${eddydir}
        cp ${rawdir}/series_index.txt ${eddydir}
        ${FSLDIR}/bin/immv ${rawdir}/Pos_Neg ${eddydir}
        mv ${rawdir}/Pos_Neg.bvals ${eddydir}
        mv ${rawdir}/Pos_Neg.bvecs ${eddydir}
        mv ${rawdir}/Pos.bv?? ${eddydir}
        mv ${rawdir}/Neg.bv?? ${eddydir}
    else
        echo ${i} diffusion data basic preprocessing results moved to appropriate directories
    fi
    
#### Run Topup (run_topup.sh with edits) ####
    if [ ! -f ${topupdir}/nodif_brain.nii.gz ]
    then
        bash ${HCPPIPEDIR}/DiffusionPreprocessing/scripts/edit_run_topup.sh ${topupdir}
    else
        echo ${i} diffusion data Topup completed
    fi

#### Run Eddy (from run_eddy.sh) #### 
    if [ ! -f ${eddydir}/eddy_unwarped_images.nii.gz ]
    then
        source ${HCPPIPEDIR}/global/scripts/log.shlib

        fsl_version_file="${FSLDIR}/etc/fslversion"
        fsl_version=`cat ${fsl_version_file}`
        echo FSL version in use is ${fsl_version}       ## for FSL 5.0.9 and above, can run GPU-enabled version of Eddy. Current script does not include GPU-enabled Eddy (should be figured out)

        ${FSLDIR}/bin/imcp ${topupdir}/nodif_brain_mask ${eddydir}/

        stdEddy="${FSLDIR}/bin/eddy"
        eddyExec="${stdEddy}"
        log_Msg "eddy executable command to use: ${eddyExec}"

        eddy_command="${eddyExec} "
        eddy_command+="--imain=${eddydir}/Pos_Neg "
        eddy_command+="--mask=${eddydir}/nodif_brain_mask "
        eddy_command+="--index=${eddydir}/series_index.txt "
        eddy_command+="--acqp=${eddydir}/acqparams.txt "
        eddy_command+="--bvecs=${eddydir}/Pos_Neg.bvecs "
        eddy_command+="--bvals=${eddydir}/Pos_Neg.bvals "
        eddy_command+="--fwhm=0 "
        eddy_command+="--topup=${topupdir}/topup_Pos_Neg_b0 "
        eddy_command+="--out=${eddydir}/eddy_unwarped_images "
        eddy_command+="--flm=quadratic "

        log_Msg "About to issue the following eddy command: "
        log_Msg "${eddy_command}"
        ${eddy_command}
        eddyReturnValue=$?
        log_Msg "Completed with return value: ${eddyReturnValue}"
        exit ${eddyReturnValue}
    else
        echo ${i} diffusion data Eddy completed
    fi





done










