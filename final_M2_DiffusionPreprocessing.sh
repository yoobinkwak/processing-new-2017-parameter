HCPPIPEDIR=/Users/yoobin_kwak/preproc_BCS_001_002/Pipelines-master
HCPPIPEDIR_Global=${HCPPIPEDIR}/global/scripts
GlobalScripts=${HCPPIPEDIR_Global} 
HCPPIPEDIR_Config=${HCPPIPEDIR}/global/config
HCPPIPEDIR_dMRI=${HCPPIPEDIR}/DiffusionPreprocessing/scripts
FSLDIR=/usr/local/fsl
FREESURFER_HOME=/Applications/freesurfer
CARET7DIR=/Users/yoobin_kwak/preproc_BCS_001_002/workbench/bin_macosx64

for i in $@
do
    Subject=${i}
    DWIName=Preprocess_Diffusion

    outdir=${Subject}/${DWIName}
    if [ ! -d ${outdir} ] ; then
        mkdir ${outdir}
    fi
    outdirT1w=${Subject}/Preprocess_Structure/PreFreesurfer/T1w/${DWIName}
    if [ ! -d ${outdirT1w} ] ; then
        mkdir ${outdirT1w}
    fi
    rawdir=${outdir}/rawdata
    if [ ! -d ${rawdir} ] ; then
        mkdir ${rawdir}
    fi
    topupdir=${outdir}/topup
    if [ ! -d ${topupdir} ] ; then
        mkdir ${topupdir}
    fi
    eddydir=${outdir}/eddy
    if [ ! -d ${eddydir} ] ; then
        mkdir ${eddydir}
    fi
    datadir=${outdir}/data
    if [ ! -d ${datadir} ] ; then
        mkdir ${datadir}
    fi
    regdir=${outdir}/reg
    if [ ! -d ${regdir} ] ; then
        mkdir ${regdir}
    fi

    echo_spacing=0.87       ## in msec
    PEdir=1         ## RL/LR phase encoding (2 if PA/AP phase encoding)
    b0dist=45       ## Minimum distance in volums between b0s considered for preprocessing (default)
    b0maxbval=50        # DEFAULT_B0_MAX_BVAL=50 from HCP 
    DegreesOfFreedom=6      ## DEFAULT_DEGREES_OF_FREEDOM=6 from HCP 
    GdCoeffs=NONE       ## file containing coefficients that describe spatial variations of the scanner gradients (need to get it from Simens)
    CombineDataFlag=2        ## If JAC resampling has been used in eddy, this value determines what to do with the output file 
                                ## 2 - include in the output all volumes uncombined (i.e. output file of eddy)
                                ## 1 - include in the output and combine only volumes where both LR/RL (or AP/PA) pairs have been acquired (Default)
                                ## 0 - As 1, but also include uncombined single volumes
    if [ ${PEdir} -eq 1 ]; then 
        basePos="RL"
        baseNeg="LR"
    elif [ ${PEdir} -eq 2 ]; then
        basePos="PA"
        baseNeg="AP"
    fi

#### Copy raw data to working directory ####
    if [ ! -f ${rawdir}/RL_BLIP_input.nii.gz ] ; then
        echo "${i} Diffusion Basic Preprocessing : Copying raw data to working directory"
        cp -r ${i}/DTI_BLIP_RL*/*.nii.gz ${rawdir}/${basePos}_BLIP_input.nii.gz       
        cp -r ${i}/DTI_BLIP_RL*/*.bval ${rawdir}/${basePos}_BLIP_input.bval 
        cp -r ${i}/DTI_BLIP_RL*/*.bvec ${rawdir}/${basePos}_BLIP_input.bvec 
        cp -r ${i}/DTI_BLIP_LR*/*.nii.gz ${rawdir}/${baseNeg}_BLIP_input.nii.gz 
        cp -r ${i}/DTI_BLIP_LR*/*.bval ${rawdir}/${baseNeg}_BLIP_input.bval 
        cp -r ${i}/DTI_BLIP_LR*/*.bvec ${rawdir}/${baseNeg}_BLIP_input.bvec
        cp -r ${i}/DTI_MB3_LR_B1000*/*.nii.gz ${rawdir}/${baseNeg}_EMPTY_B1000_input.nii.gz
        cp -r ${i}/DTI_MB3_LR_B1000*/*.bval ${rawdir}/${baseNeg}_EMPTY_B1000_input.bval
        cp -r ${i}/DTI_MB3_LR_B1000*/*.bvec ${rawdir}/${baseNeg}_EMPTY_B1000_input.bvec
        cp -r ${i}/DTI_MB3_LR_B2000*/*.nii.gz ${rawdir}/${baseNeg}_EMPTY_B2000_input.nii.gz
        cp -r ${i}/DTI_MB3_LR_B2000*/*.bval ${rawdir}/${baseNeg}_EMPTY_B2000_input.bval
        cp -r ${i}/DTI_MB3_LR_B2000*/*.bvec ${rawdir}/${baseNeg}_EMPTY_B2000_input.bvec
        cp -r ${i}/DTI_MB3_LR_B3000*/*.nii.gz ${rawdir}/${baseNeg}_EMPTY_B3000_input.nii.gz
        cp -r ${i}/DTI_MB3_LR_B3000*/*.bval ${rawdir}/${baseNeg}_EMPTY_B3000_input.bval
        cp -r ${i}/DTI_MB3_LR_B3000*/*.bvec ${rawdir}/${baseNeg}_EMPTY_B3000_input.bvec
    fi
    if [ -e ${rawdir}/RL_BLIP_input.nii.gz ] ; then
        raw_inputs=15
        check_raw_inputs=`ls ${rawdir}/*input.* | wc -l`
        if [ $raw_inputs -eq $check_raw_inputs ] ; then
            echo "${i} Diffusion Basic Preprocessing : Copied raw data to working directory"
        else
            echo "${i} Diffusion Basic Preprocessing : ERROR Copying raw data to working directorty"
        fi
    fi

#### Create two files for each phase encoding direction, that for each series contain the number of corresponding volumes and the number of actual volumes ####
    # *_SeriesCorrespVolNum.txt contains as many rows as non-EMPTY series. The entry M in row J indicates that volumes 0-M from RLseries J has corresponding LR pairs. This file is used in basic_preproc to generate topup/eddy indices and extract corresponding b0s for topup 
    # *_SeriesVolNum.txt has as many rows as maximum series pairs (even unmatched pairs). The entry M N in row J indicates that the RLSeries J has its 0-M volumes corresponding to LRSeries J and RLJ has N  volumes in total. This file is used in eddy_combine
#    MissingFileFlag="EMPTY"
#
#    Pos_count=1
#    for Image in ${rawdir}/${basePos}_BLIP.nii.gz ; do
#        if [[ ${Image} =~ ^.*EMPTY.*$  ]] ; then
#            Image=EMPTY
#        fi
#        if [ ${Image} = ${MissingFileFlag} ] ; then
#            PosVols[${Pos_count}]=0
#        else
#            PosVols[${Pos_count}]=`${FSLDIR}/bin/fslval ${Image} dim4`
#        fi
#        Pos_count=$((${Pos_count} + 1))
#    done
#
#    Neg_count=1
#    for Image in ${rawdir}/*${baseNeg}*.nii.gz ; do
#        if [[ ${Image} =~ ^.*EMPTY.*$  ]] ; then
#            Image=EMPTY
#        fi
#        if [ ${Image} = ${MissingFileFlag} ] ; then
#            NegVols[${Neg_count}]=0
#        else
#            NegVols[${Neg_count}]=`${FSLDIR}/bin/fslval ${Image} dim4`
#        fi
#        Neg_count=$((${Neg_count} + 1))
#    done
#
#    Paired_flag=0
#    for (( j=1; j<${Pos_count}; j++ )) ; do
#        CorrVols=`min ${NegVols[${j}]} ${PosVols[${j}]}`       #### NEED TO FIGURE THIS PART OUT
#        echo ${CorrVols} ${PosVols[${j}]} >> ${eddydir}/Pos_SeriesVolNum.txt
#        if [ ${PosVols[${j}]} -ne 0 ] ;then
#            echo ${CorrVols} >> ${rawdir}/${basePos}_SeriesCorrespVolNum.txt
#            if [ ${CorrVols} -ne 0 ] ; then
#                Paired_flag=1
#            fi
#        fi
#    done
#    for (( j=1; j<${Neg_count}; j++ )) ; do
#        CorrVols=`min ${NegVols[${j}]} ${PosVols[${j}]}`       #### NEED TO FIGURE THIS PART OUT
#        echo ${CorrVols} ${NegVols[${j}]} >> ${eddydir}/Neg_SeriesVolNum.txt
#        if [ ${NegVols[${j}]} -ne 0 ] ;then
#            echo ${CorrVols} >> ${rawdir}/${baseNeg}_SeriesCorrespVolNum.txt
#        fi
#    done
    
#### Basic Preprocessing (basic_preproc.sh) ####
    grappa_factor=3
    
    ##### Compute Total_readout in secs with up to 6 decimal places ####
    dimP=`${FSLDIR}/bin/fslval ${rawdir}/RL_BLIP_input.nii.gz dim1`
    nPEsteps=$(($dimP - 1))
    ro_time=`echo "${echo_spacing} * (${nPEsteps}/${grappa_factor})" | bc -l`
    ro_time=`echo "scale=6; ${ro_time} / 1000" | bc -l`
    echo "${i} Diffusion Basic Preprocessing : Total readout time is $ro_time secs"

    #### Intensity Normalisation across Series ####
    if [ ! -f ${rawdir}/LR_BLIP_input_orig.nii.gz ] ; then
        echo "${i} Diffusion Basic Preprocessing : Rescaling series to ensure consistency across baseline intensities"
        entry_cnt=0
        for entry in ${rawdir}/${basePos}*_input.nii* ${rawdir}/${baseNeg}*_input.nii* ; do
            basename=`imglob ${entry}`
            echo "      Processing $basename"
        
            echo "      $basename Create mean file" 
            if [ ! -e ${basename}_mean.nii.gz ] ; then
                ${FSLDIR}/bin/fslmaths ${entry} -Xmean -Ymean -Zmean ${basename}_mean
            fi

            echo "      ${basename} Get Posbvals from ${basename}.bval"
            Posbvals=`cat ${basename}.bval`
            echo "          Posbvals: ${Posbvals}"

            echo "      ${basename} extract all b0s"
            mcnt=0
            for p in ${Posbvals} ; do
                echo "          Posbvals p: ${p}"
                cnt=`$FSLDIR/bin/zeropad $mcnt 4`
                echo "          cnt: ${cnt}"
                if [ $p -lt ${b0maxbval} ]; then
                    $FSLDIR/bin/fslroi ${basename}_mean ${basename}_b0_${cnt} ${mcnt} 1
                fi
                mcnt=$((${mcnt} + 1))
            done

            ${FSLDIR}/bin/fslmerge -t ${basename}_mean `echo ${basename}_b0_????.nii*`
            ${FSLDIR}/bin/fslmaths ${basename}_mean -Tmean ${basename}_mean         #This is the mean baseline b0 intensity for the series
            #${FSLDIR}/bin/imrm ${basename}_b0_????

            if [ ${entry_cnt} -eq 0 ]; then      #Do not rescale the first series
                rescale=`fslmeants -i ${basename}_mean`
            else
                scaleS=`fslmeants -i ${basename}_mean`
                ${FSLDIR}/bin/fslmaths ${basename} -mul ${rescale} -div ${scaleS} ${basename}_new
                ${FSLDIR}/bin/immv ${basename} ${basename}_orig     # original data (i.e., not rescaled); HCP removes this data (imrm ${basename})
                ${FSLDIR}/bin/immv ${basename}_new ${basename}      # "replace" the original dataseries with the rescaled one
            fi
            entry_cnt=$((${entry_cnt} + 1))
            #${FSLDIR}/bin/imrm ${basename}_mean
        done
    else
        echo "${i} Diffusion Basic Preprocessing : b0 Intesnsity Normalized"
    fi

    #### b0 extraction and Creation of Index files for topup/eddy ####
    correct_Pos_b0Extract=1
    correct_Neg_b0Extract=4
    correct_b0Extract=$(($correct_Pos_b0Extract + $correct_Neg_b0Extract))
    correct_index=131
    Pos_index=7

    scount=1
    indcount=0
#    scount2=1

    if [ ! -f ${rawdir}/Pos_b0_0000.nii.gz ] ; then
        echo "${i} Diffusion Basic Preprocessing : Extracting b0s from PE_Positive volumes,  writing acqparams.txt, and creating series_index.txt (need to work on index file)"
        declare -i sesdimt
#       tmp_indx=1
#       while read line ; do  #Read SeriesCorrespVolNum.txt file
#           PCorVolNum[${tmp_indx}]=`echo $line | awk {'print $1'}`
#           tmp_indx=$((${tmp_indx}+1))
#       done < ${rawdir}/${basePos}_SeriesCorrespVolNum.txt
       # for entry in ${rawdir}/RL_BLIP_input.nii.gz ; do
        for entry in ${rawdir}/${basePos}*_input.nii* ; do
            basename=`imglob ${entry}`
            Posbvals=`cat ${basename}.bval`
            count=0     #Within series counter
            count3=$((${b0dist} + 1))
            for P in ${Posbvals} ; do
#               if [ $count -ge ${PCorVolNum[${scount2}]} ]; then
#                   tmp_ind=${indcount}
#                   if [ $[tmp_ind] -eq 0 ]; then
#                   tmp_ind=$((${indcount}+1))
#                   fi
#                   echo ${tmp_ind} >>${rawdir}/index.txt
#               else
                if [ ! -f ${rawdir}/acqparams.txt ] ; then
                    if [ $P -lt ${b0maxbval} ] && [ ${count3} -gt ${b0dist} ]; then         #Consider a b=0 a volume that has a bvalue<50 and is at least 50 volumes away from the previous
                        cnt=`$FSLDIR/bin/zeropad $indcount 4`
                        echo "Extracting Pos Volume $count from ${entry} as a b=0. Measured b=$P" >>${rawdir}/extractedb0.txt
                        $FSLDIR/bin/fslroi ${entry} ${rawdir}/Pos_b0_${cnt} ${count} 1
                        if [ ${PEdir} -eq 1 ]; then    #RL/LR phase encoding
                            echo 1 0 0 ${ro_time} >> ${rawdir}/acqparams.txt
                        fi
                        indcount=$((${indcount} + 1))
                        count3=0
                    fi
#                   echo ${indcount} >>${rawdir}/index.txt
                    count3=$((${count3} + 1))
                fi
                count=$((${count} + 1))
            done
            check_acqparams=`cat ${rawdir}/acqparams.txt | wc -l` 
            check_extractedb0=`cat ${rawdir}/extractedb0.txt |wc -l`
            if [ $check_acqparams -eq $correct_Pos_b0Extract ] && [ $check_extractedb0 -eq $correct_Pos_b0Extract ] ; then
                echo "      Extracted b0 and Wrote acqparams.txt from PE_Positive volume"
            else
                echo "      ERROR Extracting b0 writing acqparams.txt from PE_Positive volumes"
            fi
            #### Create series file ####
            if [ ! -f ${rawdir}/series_index.txt ] ; then
                sesdimt=`${FSLDIR}/bin/fslval ${entry} dim4` #Number of datapoints per Pos series
                for (( j=0; j<${sesdimt}; j++ )) ; do
                    echo ${scount} >> ${rawdir}/series_index.txt
                done
                scount=$((${scount} + 1))
#                scount2=$((${scount2} + 1))
            fi
            check_index=`cat ${rawdir}/series_index.txt | wc -l`
            if [ $check_index -eq $Pos_index ] ; then
                echo "      Created series_index.txt from PE_Positive volumes"
            else
                echo "      ERROR EXTRACTING series_index.txt from PE_Positive volumes"
            fi
        done
    fi
    
    if [ ! -f ${rawdir}/Neg_b0_0001.nii.gz ] ; then
        echo "${i} Diffusion Basic Preprocessing : Extracting b0s from PE_Negative volumes,  writing acqparams.txt, and creating series_index.txt (need to work on index file)"
        tmp_indx=1
        while read line ; do  #Read SeriesCorrespVolNum.txt file
            NCorVolNum[${tmp_indx}]=`echo $line | awk {'print $1'}`
            tmp_indx=$((${tmp_indx}+1))
        done < ${rawdir}/${baseNeg}_SeriesCorrespVolNum.txt
        Poscount=${indcount}
        for entry in ${rawdir}/${baseNeg}*_input.nii* ; do
            basename=`imglob ${entry}`
            Negbvals=`cat ${basename}.bval`
            count=0
            count3=$((${b0dist} + 1))
            for N in ${Negbvals} ; do
#                if [ $count -ge ${NCorVolNum[${scount2}]} ]; then
#                    tmp_ind=${indcount}
#                    if [ $[tmp_ind] -eq 0 ]; then
#                        tmp_ind=$((${indcount}+1))
#                    fi
#                    echo $((${tmp_ind} + ${Poscount})) >>${rawdir}/index.txt
#                else
                correct_inputs=5
                check_acqparams=`cat ${rawdir}/acqparams.txt | wc -l` 
                check_extractedb0=`cat ${rawdir}/extractedb0.txt |wc -l`
                if [ $check_acqparams -eq $correct_inputs ] && [ $check_extractedb0 -eq $correct_inputs ] ; then
                    echo "${i} Diffusion Basic Preprocessing : b0 extraction and acqparam.txt completed : to be used for topup/eddy"
                elif [ $check_acqparams -gt $correct_inputs ] && [ $check_extractedb0 -gt $correct_inputs ] ;then
                    echo "      ERROR EXTRACTING acqparams.txt and extractedb0.txt"
                else
                    if [ $N -lt ${b0maxbval} ] && [ ${count3} -gt ${b0dist} ]; then
                         cnt=`$FSLDIR/bin/zeropad $indcount 4`
                         echo "Extracting Neg Volume $count from ${entry} as a b=0. Measured b=$N" >>${rawdir}/extractedb0.txt
                         $FSLDIR/bin/fslroi ${entry} ${rawdir}/Neg_b0_${cnt} ${count} 1
                         if [ ${PEdir} -eq 1 ]; then    #RL/LR phase encoding
                             echo -1 0 0 ${ro_time} >> ${rawdir}/acqparams.txt
                         fi
                         indcount=$((${indcount} + 1))
                         count3=0
                     fi
#                     echo $((${indcount} + ${Poscount})) >>${rawdir}/index.txt
                     count3=$((${count3} + 1))
                 fi
                 count=$((${count} + 1))
             done
            
             #### Create series file ####
             check_index=`cat ${rawdir}/series_index.txt | wc -l`
             if [ $check_index -eq $correct_index ] ; then
                 echo "${i} Diffusion Basic Preprocessing : series_index.txt completed : to be used for eddy"
             elif [ $check_index -gt $correct_index ] ; then
                 echo "      ERROR EXTRACTING series_index.txt"
             else 
                 sesdimt=`${FSLDIR}/bin/fslval ${entry} dim4`
                 for (( j=0; j<${sesdimt}; j++ )) ; do
                     echo ${scount} >> ${rawdir}/series_index.txt #Create series file
                 done
                 scount=$((${scount} + 1))
                 scount2=$((${scount2} + 1)) 
             fi
         done
     fi

     if [ -e ${rawdir}/Pos_b0_0000.nii.gz ] && [ -e ${rawdir}/Neg_b0_0001.nii.gz ] && [ -e ${rawdir}/acqparams.txt ] && [ -e ${rawdir}/series_index.txt ] ; then
         Pos_b0Extracted=`ls -l ${rawdir}/Pos_b0_000*.nii.gz | wc -l`
         Neg_b0Extracted=`ls -l ${rawdir}/Neg_b0_000*.nii.gz | wc -l`
         check_acqparams=`cat ${rawdir}/acqparams.txt | wc -l` 
         check_index=`cat ${rawdir}/series_index.txt | wc -l`
         check_extractedb0=`cat ${rawdir}/extractedb0.txt |wc -l`
         if [ $correct_Pos_b0Extract -eq $Pos_b0Extracted ] &&  [ $correct_Neg_b0Extract -eq $Neg_b0Extracted ] && [ $check_acqparams -eq $correct_b0Extract ] && [ $check_index -eq $correct_index ] ; then
             echo "${i} Diffusion Basic Preprocessing : Extracted b0s, wrote acqparams.txt, and created series_index.txt from PE_Positive and PE_Negative volumes (need to work on index file)"
         else
             echo "${i} Diffusion Basic Preprocessing : ERROR Extracting b0s, writing acqparams.txt, and creating series_index.txt from PE_Positive and PE_Negative volumes (need to work on index file)"
         fi
     fi

     ##### Merging Files and correct number of slices ####
     if  [ ! -f ${rawdir}/Pos.nii.gz ] && [ ! -f ${rawdir}/Neg.nii.gz ] ; then
        echo "${i} Difussion Preprocessing : Merging Pos and Neg Images"
        ${FSLDIR}/bin/fslmerge -t ${rawdir}/Pos_b0 `${FSLDIR}/bin/imglob ${rawdir}/Pos_b0_????.*`
        ${FSLDIR}/bin/fslmerge -t ${rawdir}/Neg_b0 `${FSLDIR}/bin/imglob ${rawdir}/Neg_b0_????.*`
#        ${FSLDIR}/bin/imrm ${rawdir}/Pos_b0_????
#        ${FSLDIR}/bin/imrm ${rawdir}/Neg_b0_????
        ${FSLDIR}/bin/fslmerge -t ${rawdir}/Pos `echo ${rawdir}/${basePos}*_input.nii*`
        ${FSLDIR}/bin/fslmerge -t ${rawdir}/Neg `echo ${rawdir}/${baseNeg}*_input.nii*` 
        paste ${rawdir}/${basePos}*.bval >${rawdir}/Pos.bval
        paste ${rawdir}/${basePos}*.bvec >${rawdir}/Pos.bvec
        paste ${rawdir}/${baseNeg}*.bval >${rawdir}/Neg.bval
        paste ${rawdir}/${baseNeg}*.bvec >${rawdir}/Neg.bvec
    fi
    output_Pos=3
    output_Neg=3
    output_Pos_Neg=4
    if  [  -e ${rawdir}/Pos.nii.gz ] && [  -e ${rawdir}/Neg.nii.gz ] ; then
        count_Pos=`ls -l ${rawdir}/Pos.* | wc -l`
        count_Neg=`ls -l ${rawdir}/Neg.* | wc -l`
        if [ $output_Pos -eq $count_Pos ] &&  [ $output_Neg -eq $count_Neg ] && [ -e ${rawdir}/Pos_b0.nii.gz ] && [ -e ${rawdir}/Neg_b0.nii.gz ] ; then
            echo "${i} Diffusion Basic Preprocessing : Pos and Neg Images Merged"
        else
            echo "${i} Diffusion Basic Preprocessing : ERROR MERGING POS AND NEG IMAGES" 
        fi
    fi

    dimz=`${FSLDIR}/bin/fslval ${rawdir}/Pos dim3`       
    if [ $((dimz%2)) -eq 0 ] ; then
        echo "${i} Diffusion Basic Preprocessing : Slice Number $dimz is Even"
    else
        echo "${i} Diffusion Basic Preprocessing : Remove one slice from data to get even number of slices"
        ${FSLDIR}/bin/fslroi ${rawdir}/Pos ${rawdir}/Posn 0 -1 0 -1 1 -1
        ${FSLDIR}/bin/fslroi ${rawdir}/Neg ${rawdir}/Negn 0 -1 0 -1 1 -1
        ${FSLDIR}/bin/fslroi ${rawdir}/Pos_b0 ${rawdir}/Pos_b0n 0 -1 0 -1 1 -1
        ${FSLDIR}/bin/fslroi ${rawdir}/Neg_b0 ${rawdir}/Neg_b0n 0 -1 0 -1 1 -1
        ${FSLDIR}/bin/imrm ${rawdir}/Pos
        ${FSLDIR}/bin/imrm ${rawdir}/Neg
        ${FSLDIR}/bin/imrm ${rawdir}/Pos_b0
        ${FSLDIR}/bin/imrm ${rawdir}/Neg_b0
        ${FSLDIR}/bin/immv ${rawdir}/Posn ${rawdir}/Pos
        ${FSLDIR}/bin/immv ${rawdir}/Negn ${rawdir}/Neg
        ${FSLDIR}/bin/immv ${rawdir}/Pos_b0n ${rawdir}/Pos_b0
        ${FSLDIR}/bin/immv ${rawdir}/Neg_b0n ${rawdir}/Neg_b0
    fi

    if [ ! -f ${rawdir}/Pos_Neg.nii.gz ] ; then
        echo "${i} Diffusion Basic Preprocessing : Perform Final Merge"
        ${FSLDIR}/bin/fslmerge -t ${rawdir}/Pos_Neg_b0 ${rawdir}/Pos_b0 ${rawdir}/Neg_b0
        ${FSLDIR}/bin/fslmerge -t ${rawdir}/Pos_Neg ${rawdir}/Pos ${rawdir}/Neg
        paste ${rawdir}/Pos.bval ${rawdir}/Neg.bval >${rawdir}/Pos_Neg.bvals
        paste ${rawdir}/Pos.bvec ${rawdir}/Neg.bvec >${rawdir}/Pos_Neg.bvecs
#        ${FSLDIR}/bin/imrm ${rawdir}/Pos
#        ${FSLDIR}/bin/imrm ${rawdir}/Neg
    fi
    if [ -e  ${rawdir}/Pos_Neg.nii.gz ] ; then
        count_Pos=`ls -l ${rawdir}/Pos.* | wc -l`
        count_Neg=`ls -l ${rawdir}/Neg.* | wc -l`
        count_Pos_Neg=`ls -l ${rawdir}/Pos_Neg* | wc -l`
        if [ $count_Pos_Neg -eq $output_Pos_Neg ] ; then
            echo "${i} Diffusion Basic Preprocessing : Final Merge Completed"
        else
            echo "${i} Diffusion Basic Preprocessing : ERROR PERFORMING FINAL MERGE"
        fi
    fi

    ##### Move files to appropriate directories ####
    if [ -f ${rawdir}/Neg.bvec ] && [ ! -f ${eddydir}/Neg.bvec ] ; then
        echo "${i} Diffusion Basic Preprocessing : Moving Files to Appropriate Directories"
        cp ${rawdir}/extractedb0.txt ${topupdir}/extractedb0.txt
        cp ${rawdir}/acqparams.txt ${topupdir}/acqparams.txt
        ${FSLDIR}/bin/imcp ${rawdir}/Pos_Neg_b0 ${topupdir}
        ${FSLDIR}/bin/imcp ${rawdir}/Pos_b0 ${topupdir}
        ${FSLDIR}/bin/imcp ${rawdir}/Neg_b0 ${topupdir}
        cp ${topupdir}/acqparams.txt ${eddydir}
#        cp ${rawdir}/index.txt ${eddydir}
        cp ${rawdir}/series_index.txt ${eddydir}
        ${FSLDIR}/bin/imcp ${rawdir}/Pos_Neg ${eddydir}
        cp ${rawdir}/Pos_Neg.bvals ${eddydir}
        cp ${rawdir}/Pos_Neg.bvecs ${eddydir}
        cp ${rawdir}/Pos.bv?? ${eddydir}
        cp ${rawdir}/Neg.bv?? ${eddydir}
    else
        echo "${i} Diffusion Basic Preprocessing : Moved Files to Appropriate Directories"
    fi
    
##### Run Topup (from run_topup.sh ####
    topup_config_file=${FSLDIR}/etc/flirtsch/b02b0.cnf

    if [ ! -f ${topupdir}/topup_Pos_Neg_b0_fieldcoef.nii.gz ] ; then
        echo "${i} Diffusion Topup Processing : Perform Topup" 
        ${FSLDIR}/bin/topup --imain=${topupdir}/Pos_Neg_b0 --datain=${topupdir}/acqparams.txt --config=${topup_config_file} --out=${topupdir}/topup_Pos_Neg_b0 -v
    else
        echo "${i} Diffusion Topup Processing : Topup Completed"
    fi
        
    dimt=`${FSLDIR}/bin/fslval ${topupdir}/Pos_b0 dim4`
    dimt=$((${dimt} + 1))
    if [ ! -f ${topupdir}/hifib0.nii.gz ] ; then
        echo "${i} Diffusion Topup Processing : Perform Applytopup" 
        ${FSLDIR}/bin/fslroi ${topupdir}/Pos_b0 ${topupdir}/Pos_b01 0 1
        ${FSLDIR}/bin/fslroi ${topupdir}/Neg_b0 ${topupdir}/Neg_b01 0 1
        ${FSLDIR}/bin/applytopup --imain=${topupdir}/Pos_b01,${topupdir}/Neg_b01 --topup=${topupdir}/topup_Pos_Neg_b0 --datain=${topupdir}/acqparams.txt --inindex=1,${dimt} --out=${topupdir}/hifib0
    else
        echo "${i} Diffusion Topup Processing : Applytopup Completed"
    fi

    if [ ! -f ${topupdir}/nodif_brain.nii.gz ] ; then
        echo "${i} Diffusion Topup Processing : Perform BET on Applytopup Output" 
        ${FSLDIR}/bin/bet ${topupdir}/hifib0 ${topupdir}/nodif_brain -m -f 0.2
    else
        echo "${i} Diffusion Topup Processing : BET on Applytopup Output Completed"
    fi

##### Run Eddy (from run_eddy.sh) #### 
    if [ ! -f ${eddydir}/eddy_unwarped_images.nii.gz ] ; then
        echo "${i} Diffusion Eddy Processing : Perform Eddy"
        fsl_version_file="${FSLDIR}/etc/fslversion"
        fsl_version=`cat ${fsl_version_file}`
        echo "          FSL version in use is ${fsl_version} (for FSL 5.0.9 and above, can run GPU-enabled version of Eddy)"          ## for FSL 5.0.9 and above, can run GPU-enabled version of Eddy. Current script does not include GPU-enabled Eddy (should be figured     out)
        if [ ! -f ${eddydir}/nodif_brain_mask.nii.gz ] ; then
            ${FSLDIR}/bin/imcp ${topupdir}/nodif_brain_mask ${eddydir}/
        fi
        ${FSLDIR}/bin/eddy --imain=${eddydir}/Pos_Neg --mask=${eddydir}/nodif_brain_mask --index=${eddydir}/series_index.txt --acqp=${eddydir}/acqparams.txt --bvecs=${eddydir}/Pos_Neg.bvecs --bvals=${eddydir}/Pos_Neg.bvals --fwhm=0 --topup=${topupdir}/topup_Pos_Neg_b0 --out=${eddydir}/eddy_unwarped_images --flm=quadratic 
    else
        echo "${i} Diffusion Eddy Processing :  Eddy Completed"
    fi
    
#### Eddy Post Processing (from eddy_postproc.sh) ####
    #### Combine Data Flag options ####
    if [ ${CombineDataFlag} -eq 2 ]; then
        echo "${i} Diffusion Post Eddy Processing : Not Combining Eddy Output (Combine Data Flag set to 2)"
        if [ ! -f ${datadir}/data.nii.gz ] ;then
            ${FSLDIR}/bin/imcp  ${eddydir}/eddy_unwarped_images ${datadir}/data
            cp ${eddydir}/Pos_Neg.bvals ${datadir}/bvals
            cp ${eddydir}/Pos_Neg.bvecs ${datadir}/bvecs
        elif [ -e ${datadir}/data.nii.gz ] && [ -e ${datadir}/bvals ] && [ -e ${datadir}/bvecs ] ; then
            echo "      Copied Eddy Outputs Copied "data" Directory"
        else
            echo "      ERROR Copying Eddy Outputs to "data" Directory"
        fi
    else
        echo "${i} Diffusion Post Eddy Processing : JAC resampling has been used. Eddy Output is now combined (NEED TO CHECK LINES 491~539)"
        PosVols=`wc ${eddydir}/Pos.bval | awk {'print $2'}`
        NegVols=`wc ${eddydir}/Neg.bval | awk {'print $2'}`    #Split Pos and Neg Volumes
        ${FSLDIR}/bin/fslroi ${eddydir}/eddy_unwarped_images ${eddydir}/eddy_unwarped_Pos 0 ${PosVols}
        ${FSLDIR}/bin/fslroi ${eddydir}/eddy_unwarped_images ${eddydir}/eddy_unwarped_Neg ${PosVols} ${NegVols}
        ${FSLDIR}/bin/eddy_combine ${eddydir}/eddy_unwarped_Pos ${eddydir}/Pos.bval ${eddydir}/Pos.bvec ${eddydir}/Pos_SeriesVolNum.txt ${eddydir}/eddy_unwarped_Neg ${eddydir}/Neg.bval ${eddydir}/Neg.bvec ${eddydir}/Neg_SeriesVolNum.txt ${datadir} ${CombineDataFlag}

        ${FSLDIR}/bin/imrm ${eddydir}/eddy_unwarped_Pos
	    ${FSLDIR}/bin/imrm ${eddydir}/eddy_unwarped_Neg
	    cp ${datadir}/bvals ${datadir}/bvals_noRot
	    cp ${datadir}/bvecs ${datadir}/bvecs_noRot
     
	    # Divide Eddy-Rotated bvecs to Pos and Neg
	    line1=`awk 'NR==1 {print; exit}' ${eddydir}/eddy_unwarped_images.eddy_rotated_bvecs`
	    line2=`awk 'NR==2 {print; exit}' ${eddydir}/eddy_unwarped_images.eddy_rotated_bvecs`
	    line3=`awk 'NR==3 {print; exit}' ${eddydir}/eddy_unwarped_images.eddy_rotated_bvecs`   
	    Posline1=""
	    Posline2=""
	    Posline3=""
	    for ((i=1; i<=$PosVols; i++)); do
	        Posline1="$Posline1 `echo $line1 | awk -v N=$i '{print $N}'`"
	        Posline2="$Posline2 `echo $line2 | awk -v N=$i '{print $N}'`"
	        Posline3="$Posline3 `echo $line3 | awk -v N=$i '{print $N}'`"
	    done
	    echo $Posline1 > ${eddydir}/Pos_rotated.bvec
	    echo $Posline2 >> ${eddydir}/Pos_rotated.bvec
	    echo $Posline3 >> ${eddydir}/Pos_rotated.bvec

	    Negline1=""
	    Negline2=""
	    Negline3=""
	    Nstart=$((PosVols + 1 ))
	    Nend=$((PosVols + NegVols))
	    for  ((i=$Nstart; i<=$Nend; i++)); do
	        Negline1="$Negline1 `echo $line1 | awk -v N=$i '{print $N}'`"
	        Negline2="$Negline2 `echo $line2 | awk -v N=$i '{print $N}'`"
	        Negline3="$Negline3 `echo $line3 | awk -v N=$i '{print $N}'`"
	    done
	    echo $Negline1 > ${eddydir}/Neg_rotated.bvec
	    echo $Negline2 >> ${eddydir}/Neg_rotated.bvec
	    echo $Negline3 >> ${eddydir}/Neg_rotated.bvec
	
	    # Average Eddy-Rotated bvecs. Get for each direction the two b matrices, average those and then eigendecompose the average b-matrix to get the new bvec and bval.
	    # Also outputs an index file (1-based) with the indices of the input (Pos/Neg) volumes that have been retained in the output
	    ${globalscriptsdir}/average_bvecs.py ${eddydir}/Pos.bval ${eddydir}/Pos_rotated.bvec ${eddydir}/Neg.bval ${eddydir}/Neg_rotated.bvec ${datadir}/avg_data ${eddydir}/Pos_SeriesVolNum.txt ${eddydir}/Neg_SeriesVolNum.txt

	    mv ${datadir}/avg_data.bval ${datadir}/bvals
	    mv ${datadir}/avg_data.bvec ${datadir}/bvecs
	    rm -f ${datadir}/avg_data.bv??
    fi

    #### Gradient nonlinearty coeeficient ####
    if [ ! $GdCoeffs = "NONE" ] ; then
        echo "${i} Diffusion Post Eddy Processing : Correcting for gradient nonlinearities (NEED TO CHECK LINES 544~557)"
        ${FSLDIR}/bin/immv ${datadir}/data ${datadir}/data_warped
        ${globalscriptsdir}/GradientDistortionUnwarp.sh --workingdir="${datadir}" --coeffs="${GdCoeffs}" --in="${datadir}/data_warped" --out="${datadir}/data" --owarp="${datadir}/fullWarp"
        echo "Computing gradient coil tensor to correct for gradient nonlinearities"
        ${FSLDIR}/bin/calc_grad_perc_dev --fullwarp=${datadir}/fullWarp -o ${datadir}/grad_dev
        ${FSLDIR}/bin/fslmerge -t ${datadir}/grad_dev ${datadir}/grad_dev_x ${datadir}/grad_dev_y ${datadir}/grad_dev_z
        ${FSLDIR}/bin/fslmaths ${datadir}/grad_dev -div 100 ${datadir}/grad_dev #Convert from % deviation to absolute
        ${FSLDIR}/bin/imrm ${datadir}/grad_dev_?
        ${FSLDIR}/bin/imrm ${datadir}/trilinear
        ${FSLDIR}/bin/imrm ${datadir}/data_warped_vol1
        #Keep the original warped data and warp fields
        mkdir -p ${datadir}/warped
        ${FSLDIR}/bin/immv ${datadir}/data_warped ${datadir}/warped
        ${FSLDIR}/bin/immv ${datadir}/fullWarp ${datadir}/warped
        ${FSLDIR}/bin/immv ${datadir}/fullWarp_abs ${datadir}/warped
    else
        echo "${i} Diffusion Post Eddy Processing : Not correcting for gradient nonlinearities"
    fi

    #Remove negative intensity values (caused by spline interpolation) from final data
    if [ ! -f ${datadir}/nodif.nii.gz ] ; then
        echo "${i} Diffusion Post Eddy Processing : Removing Negative Intensity Values from Final Eddy Data"
        ${FSLDIR}/bin/fslmaths ${datadir}/data -thr 0 ${datadir}/data
        ${FSLDIR}/bin/bet ${datadir}/data ${datadir}/nodif_brain -m -f 0.1
        $FSLDIR/bin/fslroi ${datadir}/data ${datadir}/nodif 0 1
    elif [ -e ${datadir}/nodif.nii.gz ] && [ -e ${datadir}/nodif_brain.nii.gz ] && [ -e ${datadir}/nodif_brain_mask.nii.gz ] ; then
        echo "${i} Diffusion Post Eddy Processing : Removed Negative Intensity Values from Final Eddy Data"
    else
        echo "${i} Diffusion Post Eddy Processing : ERROR Removing Negative Intensity Values from Final Eddy Data"
    fi

#### Diffusion to Structural Registration (from DiffusionToStructural.sh) ####
    T1wOutputDirectory=${outdirT1w}
    T1wFolder=${i}/Preprocess_Structure/PreFreesurfer/T1w
    T1wImage=${T1wFolder}/T1w_acpc_dc
    T1wRestoreImage=${T1wFolder}/T1w_acpc_dc_restore
    T1wRestoreImageBrain=${T1wFolder}/T1w_acpc_dc_restore_brain
    T1wBrainImage=$T1wRestoreImageBrain
    BiasField=${T1wFolder}/BiasField_acpc_dc
    InputBrainMask=${T1wFolder}/brainmask_fs       ####from PostFS (FreeSurfer2CaretConvertAndRegisterNonlinear.sh), processing incorporated here till PostFS script is completed
    RegOutput=${outdir}/reg/Scout2T1w
    QAImage=${outdir}/reg/T1wMulEPI
    DiffRes=`${FSLDIR}/bin/fslval ${outdir}/data/data pixdim1`
    DiffRes=`printf %0.2f ${DiffRes}` 
    T1wBrainImageFile=`basename $T1wBrainImage`
    regimg="nodif"

    if [ ! -f ${regdir}/T1w_acpc_dc_restore_brain.nii.gz ] ; then
        ${FSLDIR}/bin/imcp $T1wBrainImage ${regdir}/$T1wBrainImageFile
    fi

    #### b0 FLIRT BBR and bbregister to T1w (lines 604~642 from epi_reg_dof) ####
    dof=6
    vepi=${datadir}/$regimg     ## --epi
    vrefhead=$T1wImage          ## --t1
    vrefbrain=${regdir}/$T1wBrainImageFile        ## --t1brain
    vout=${regdir}/"$regimg"2T1w_initII        ##--out
    if [ ! -f ${regdir}/"$regimg"2T1w_initII.nii.gz ] ; then
        echo "${i} Diffusion to Structural Registration : Running b0 FLIRT BBR to T1w (epi_reg_dof)"
        # create the WM segmentation
        if [ ! -f ${vout}_fast_wmseg.nii.gz ] ; then
            echo "      Running FAST segmentation"
            $FSLDIR/bin/fast -o ${vout}_fast ${vrefbrain}
            $FSLDIR/bin/fslmaths ${vout}_fast_pve_2 -thr 0.5 -bin ${vout}_fast_wmseg
        else 
            echo "      Ran FAST segmentation"
        fi
        # make a WM edge map for visualisation (good to overlay in FSLView)
        if [ ! -f ${vout}_fast_wmedge.nii.gz ] ; then
            echo "      Create WM edge for visualisation"
            $FSLDIR/bin/fslmaths ${vout}_fast_wmseg -edge -bin -mas ${vout}_fast_wmseg ${vout}_fast_wmedge
        else
            echo "      WM edge for visualisation created"
        fi
        # do a standard flirt pre-alignment
        if [ ! -f ${vout}_init.mat ] ; then
            echo "      Running FLIRT pre-alignment"
            $FSLDIR/bin/flirt -ref ${vrefbrain} -in ${vepi} -dof ${dof} -omat ${vout}_init.mat
        else
            echo "      Ran FLIRT pre-alignment"
        fi
        # run the bbr
        if [ ! -f ${vout}.mat ] ; then
            echo "      Running BBR"
            $FSLDIR/bin/flirt -ref ${vrefhead} -in ${vepi} -dof ${dof} -cost bbr -wmseg ${vout}_fast_wmseg -init ${vout}_init.mat -omat ${vout}.mat -out ${vout} -schedule ${FSLDIR}/etc/flirtsch/bbr.sch
            $FSLDIR/bin/applywarp -i ${vepi} -r ${vrefhead} -o ${vout} --premat=${vout}.mat --interp=spline
        else
            echo "      Ran BBR"
        fi
    else
        echo "${i} Diffusion to Structural Registration : b0 FLIRT BBR to T1w (epi_reg_dof) Completed"
    fi

    if [ ! -f ${regdir}/"$regimg"2T1w_restore_initII.nii.gz ] ; then
        echo "      Process FLIRT BBR outputs for fine-tuning with bbregister"
        ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${datadir}/"$regimg" -r $T1wImage --premat=${regdir}/"$regimg"2T1w_initII_init.mat -o ${regdir}/"$regimg"2T1w_init.nii.gz
        ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${datadir}/"$regimg" -r $T1wImage --premat=${regdir}/"$regimg"2T1w_initII.mat -o ${regdir}/"$regimg"2T1w_initII.nii.gz
        ${FSLDIR}/bin/fslmaths ${regdir}/"$regimg"2T1w_initII.nii.gz -div $BiasField ${regdir}/"$regimg"2T1w_restore_initII.nii.gz
    else
        echo "      FLIRT BBR outputs Processed for fine-tuning with bbregister"
    fi

    FreeSurferSubjectFolder=${i}/Preprocess_Structure/Freesurfer
    FreeSurferSubjectID=${Subject}
    SUBJECTS_DIR=$FreeSurferSubjectFolder
    export SUBJECTS_DIR
    if [ ! -f ${regdir}/"$regimg"2T1w.nii.gz ] ; then
        echo "${i} Diffusion to Structural Registration : Fine-tuning b0 to T1w Registration with bbregister"
        ${FREESURFER_HOME}/bin/bbregister --s $FreeSurferSubjectID --mov ${regdir}/"$regimg"2T1w_restore_initII.nii.gz --surf white.deformed --init-reg $FreeSurferSubjectFolder/$FreeSurferSubjectID/mri/transforms/eye.dat --bold --reg ${regdir}/EPItoT1w.dat --o ${regdir}/"$regimg"2T1w.nii.gz
        if [ ! -f ${regdir}/diff2str_fs.mat ] ; then
            echo "${i} Diffusion to Structural Registration : Running tkregister2"
            ${FREESURFER_HOME}/bin/tkregister2 --noedit --reg ${regdir}/EPItoT1w.dat --mov ${regdir}/"$regimg"2T1w_restore_initII.nii.gz --targ "$T1wImage".nii.gz --fslregout ${regdir}/diff2str_fs.mat
        fi
    elif [  -e ${regdir}/"$regimg"2T1w.nii.gz ] && [ ! -f ${regdir}/diff2str_fs.mat ] ; then
        echo "${i} Diffusion to Structural Registration : Running tkregister2"
        ${FREESURFER_HOME}/bin/tkregister2 --noedit --reg ${regdir}/EPItoT1w.dat --mov ${regdir}/"$regimg"2T1w_restore_initII.nii.gz --targ "$T1wImage".nii.gz --fslregout ${regdir}/diff2str_fs.mat
    elif [  -e ${regdir}/"$regimg"2T1w.nii.gz ] && [ -e ${regdir}/diff2str_fs.mat ] ; then
        echo "${i} Diffusion to Structural Registration : b0 to T1w Registration Fine-tuned with bbregister"
    else
        echo "${i} Diffusion to Structural Registration : ERROR Fine-tuning b0 to T1w Registration with bbregister"
    fi

    if [ ! -f ${regdir}/diff2str.mat ] ; then
        ${FSLDIR}/bin/convert_xfm -omat ${regdir}/diff2str.mat -concat ${regdir}/diff2str_fs.mat ${regdir}/"$regimg"2T1w_initII.mat
    fi
    if [ ! -f ${regdir}/str2diff.mat ] ; then
        ${FSLDIR}/bin/convert_xfm -omat ${regdir}/str2diff.mat -inverse ${regdir}/diff2str.mat
    fi

    if [ ! -f ${regdir}/"$regimg"2T1w.nii.gz ] || [ ! -f ${regdir}/"$regimg"2T1w_restore.nii.gz ] ; then
        echo "${i} Diffusion to Structural Registration : Warping with tkregister2 output"
        ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${datadir}/$regimg -r "$T1wImage".nii.gz --premat=${regdir}/diff2str.mat -o ${regdir}/"$regimg"2T1w
        ${FSLDIR}/bin/fslmaths ${regdir}/"$regimg"2T1w -div $BiasField ${regdir}/"$regimg"2T1w_restore
    elif [ -e ${regdir}/"$regimg"2T1w.nii.gz ] && [ -e ${regdir}/"$regimg"2T1w_restore.nii.gz ] ; then
        echo "${i} Diffusion to Structural Registration : Warped with tkregister2 output"
    else
        echo "${i} Diffusion to Structural Registration : ERROR Warping with tkregister2 output"
    fi
    #Are the next two scripts needed?
    if [ ! -f ${outdir}/reg/Scout2T1w.nii.gz ] ; then
        ${FSLDIR}/bin/imcp ${regdir}/"$regimg"2T1w_restore $RegOutput
    fi
    if [ ! -f ${outdir}/reg/T1wMulEPI_nodif.nii.gz ] ; then
        ${FSLDIR}/bin/fslmaths "$T1wRestoreImage".nii.gz -mul ${regdir}/"$regimg"2T1w_restore.nii.gz -sqrt "$QAImage"_"$regimg".nii.gz
    fi

    ##### Generate 2.3mm structural space for resampling the diffusion data into (1.25mm for HCP) ####     
    if [ ! -f ${T1wFolder}/T1w_acpc_dc_restore_${DiffRes}.nii.gz ] ; then
        echo "${i} Diffusion to Structural Registration : Generating ${DiffRes} native structural volume space"
        ${FSLDIR}/bin/flirt -interp spline -in $T1wRestoreImage -ref $T1wRestoreImage -applyisoxfm ${DiffRes} -out "$T1wRestoreImage"_${DiffRes}
        ${FSLDIR}/bin/applywarp --rel --interp=spline -i $T1wRestoreImage -r "$T1wRestoreImage"_${DiffRes} -o "$T1wRestoreImage"_${DiffRes}
    elif [ -e ${T1wFolder}/T1w_acpc_dc_restore_${DiffRes}.nii.gz ] ; then
        echo "${i} Diffusion to Structural Registration : ${DiffRes} Native Structural Volume Space Generated"
    else
        echo "${i} Diffusion to Structural Registration :ERROR Generating ${DiffRes} native structural volume space"
    fi
    
    ##### Generate 2.3mm mask in structural space (1.25mm for HCP) ####
    ## This step requires PostFreesurfer Output
    if [ ! -f $InputBrainMask.nii.gz ] || [ ! -f ${i}/Preprocess_Structure/PreFreesurfer/MNINonLinear/brainmask_fs.nii.gz ] ; then
        echo "${i} Diffusion to Structural Registration : Generating Brain Mask (Post Freesurfer Output)"
        AtlasSpaceFolder=${i}/Preprocess_Structure/PreFreesurfer/MNINonLinear
        AtlasSpaceT1wImage=$AtlasSpaceFolder/T1w_restore
        AtlasTransform=$AtlasSpaceFolder/xfms/acpc_dc2standard
        FreeSurferLabels=$HCPPIPEDIR/global/config/FreeSurferCorticalLabelTableLut.txt       ## NEED TO CONFIRM WHETHER THIS IS THE RIGHT FILE
        T1wImageBrainMask=brainmask_fs
        #### Convert FreeSurfer Volumes ####
        for Image in wmparc aparc.a2009s+aseg aparc+aseg ; do
            if [ ! -f $AtlasSpaceFolder/$Image.nii.gz ] ; then
                if [ -e $FreeSurferSubjectFolder/$FreeSurferSubjectID/mri/$Image.mgz ] ; then
                    echo "      Converting FreeSurfer $Image to nifti file"
                    mri_convert -rt nearest -rl $T1wImage.nii.gz $FreeSurferSubjectFolder/$FreeSurferSubjectID/mri/$Image.mgz $T1wFolder/"$Image"_1mm.nii.gz
                    applywarp --rel --interp=nn -i $T1wFolder/"$Image"_1mm.nii.gz -r $AtlasSpaceT1wImage --premat=$FSLDIR/etc/flirtsch/ident.mat -o $T1wFolder/$Image.nii.gz
                    applywarp --rel --interp=nn -i $T1wFolder/"$Image"_1mm.nii.gz -r $AtlasSpaceT1wImage -w $AtlasTransform -o $AtlasSpaceFolder/$Image.nii.gz
                    ${CARET7DIR}/wb_command -volume-label-import $T1wFolder/$Image.nii.gz $FreeSurferLabels $T1wFolder/$Image.nii.gz -drop-unused-labels
                    ${CARET7DIR}/wb_command -volume-label-import $AtlasSpaceFolder/$Image.nii.gz $FreeSurferLabels $AtlasSpaceFolder/$Image.nii.gz -drop-unused-labels
                else
                    echo "      ERRROR Converting FreesSurfer $Image to nifti file (FreeSurfer Incomplete)"
                fi
            else
                echo "      Converted FreesSurfer $Image to nifti file"
            fi
        done
        #### Create FreeSurfer Brain Mask ####
        if [ ! -f $AtlasSpaceFolder/$T1wImageBrainMask.nii.gz ] || [ ! -f $T1wFolder/$T1wImageBrainMask.nii.gz ] ; then
            echo "      Creating FreeSurfer Brain Mask"
            ${FSLDIR}/bin/fslmaths $T1wFolder/wmparc_1mm.nii.gz -bin -dilD -dilD -dilD -ero -ero $T1wFolder/"$T1wImageBrainMask"_1mm.nii.gz
            ${CARET7DIR}/wb_command -volume-fill-holes $T1wFolder/"$T1wImageBrainMask"_1mm.nii.gz $T1wFolder/"$T1wImageBrainMask"_1mm.nii.gz
            ${FSLDIR}/bin/fslmaths $T1wFolder/"$T1wImageBrainMask"_1mm.nii.gz -bin $T1wFolder/"$T1wImageBrainMask"_1mm.nii.gz
            ${FSLDIR}/bin/applywarp --rel --interp=nn -i $T1wFolder/"$T1wImageBrainMask"_1mm.nii.gz -r $AtlasSpaceT1wImage --premat=$FSLDIR/etc/flirtsch/ident.mat -o $T1wFolder/$T1wImageBrainMask.nii.gz
            ${FSLDIR}/bin/applywarp --rel --interp=nn -i $T1wFolder/"$T1wImageBrainMask"_1mm.nii.gz -r $AtlasSpaceT1wImage -w $AtlasTransform -o $AtlasSpaceFolder/$T1wImageBrainMask.nii.gz
        elif [  -e $AtlasSpaceFolder/$T1wImageBrainMask.nii.gz ] && [ -e $T1wFolder/$T1wImageBrainMask.nii.gz ] ; then
            echo "      Created FreeSurfer Brain Mask"
        else
            echo "      ERROR Creating FreeSurfer Brain Mask"
        fi
    elif [ -e $InputBrainMask.nii.gz ] && [ -e ${i}/Preprocess_Structure/PreFreesurfer/MNINonLinear/brainmask_fs.nii.gz ] ; then
        echo "${i} Diffusion to Structural Registration : Brain Mask (Post Freesurfer Output) Generated"
    else
        echo "${i} Diffusion to Structural Registration : ERROR Generating Brain Mask (Post Freesurfer Output)"
    fi

    if [ -e $InputBrainMask.nii.gz ] && [ ! -f $T1wOutputDirectory/nodif_brain_mask.nii.gz ] ; then
        echo "${i} Diffusion to Structural Registration : Generating ${DiffRes} native structural volume space mask"
        ${FSLDIR}/bin/flirt -interp nearestneighbour -in $InputBrainMask -ref $InputBrainMask -applyisoxfm ${DiffRes} -out $T1wOutputDirectory/nodif_brain_mask
        ${FSLDIR}/bin/fslmaths $T1wOutputDirectory/nodif_brain_mask -kernel 3D -dilM $T1wOutputDirectory/nodif_brain_mask
    elif [ -e $InputBrainMask.nii.gz ] && [ -e $T1wOutputDirectory/nodif_brain_mask.nii.gz ] ; then
        echo "${i} Diffusion to Structural Registration : ${DiffRes} Native Structural Volume Space Mask Generated"
    else
        echo "${i} Diffusion to Structural Registration : ERROR Generating ${DiffRes} native structural volume space mask"
    fi
    
    #### Dilated mask for masking the final data and grad_dev ####
    if [ ! -f $T1wOutputDirectory/nodif_brain_mask_temp.nii.gz ] ; then
        echo "${i} Diffusion to Structural Registration : Dialating ${DiffRes} native structural volume space mask"
        DilationsNum=6 
        ${FSLDIR}/bin/imcp $T1wOutputDirectory/nodif_brain_mask $T1wOutputDirectory/nodif_brain_mask_temp
        for (( j=0; j<${DilationsNum}; j++ )) ; do
            ${FSLDIR}/bin/fslmaths $T1wOutputDirectory/nodif_brain_mask_temp -kernel 3D -dilM $T1wOutputDirectory/nodif_brain_mask_temp
        done
    elif [ -e $T1wOutputDirectory/nodif_brain_mask_temp.nii.gz ] ; then
        echo "${i} Diffusion to Structural Registration : ${DiffRes} Native Structural Volume Space Mask Dialated"
    else
        echo "${i} Diffusion to Structural Registration : ERROR Dialating ${DiffRes} native structural volume space mask"
    fi

    #### Rotate bvecs from diffusion to structural space ####
    if [ ! -f $T1wOutputDirectory/bvecs ] || [ ! -f $T1wOutputDirectory/bvals ] ; then
        echo "${i} Diffusion to Structural Registration : Rotate bvecs from diffusion to structural space"
        ${GlobalScripts}/Rotate_bvecs.sh ${datadir}/bvecs ${regdir}/diff2str.mat $T1wOutputDirectory/bvecs
        cp ${datadir}/bvals $T1wOutputDirectory/bvals
    elif [ -e $T1wOutputDirectory/bvecs ] && [ -e $T1wOutputDirectory/bvals ] ; then
        echo "${i} Diffusion to Structural Registration : bvecs Roated from Dffusion to Structural space"
    fi

    #### Register diffusion data to T1w space. Account for gradient nonlinearities if requested ####
    if [ ! -f $T1wOutputDirectory/data.nii.gz ] ; then
        echo "${i} Diffusion to Structural Registration : Register diffusion data to T1w space"
        #### Determine whether Gradient Nonlinearity Distortion coefficients are supplied ####
        GdFlag=0
        if [ ! ${GdCoeffs} = "NONE" ] ; then
            GdFlag=1
        fi

        if [ ${GdFlag} -eq 1 ]; then
            echo "      Correcting Diffusion data for gradient nonlinearities and registering to structural space"
            ${FSLDIR}/bin/convertwarp --rel --relout --warp1=${datadir}/warped/fullWarp --postmat=${regdir}/diff2str.mat --ref="$T1wRestoreImage"_${DiffRes} --out=${regdir}/grad_unwarp_diff2str
            ${FSLDIR}/bin/applywarp --rel -i ${datadir}/warped/data_warped -r "$T1wRestoreImage"_${DiffRes} -w ${regdir}/grad_unwarp_diff2str --interp=spline -o $T1wOutputDirectory/data
            #Now register the grad_dev tensor 
            ${FSLDIR}/bin/vecreg -i ${datadir}/grad_dev -o $T1wOutputDirectory/grad_dev -r "$T1wRestoreImage"_${DiffRes} -t ${regdir}/diff2str.mat --interp=spline
            ${FSLDIR}/bin/fslmaths $T1wOutputDirectory/grad_dev -mas $T1wOutputDirectory/nodif_brain_mask_temp $T1wOutputDirectory/grad_dev  #Mask-out values outside the brain 
        else
            echo "      Registering diffusion data to T1w space without considering gradient nonlinearities"
            ${FSLDIR}/bin/flirt -in ${datadir}/data -ref "$T1wRestoreImage"_${DiffRes} -applyxfm -init ${regdir}/diff2str.mat -interp spline -out $T1wOutputDirectory/data
        fi
    else
        echo "${i} Diffusion to Structural Registration : Diffusion Data Registered to T1w Space with GdCorrection ${GdFlag} (1: On, 0: Off)"
    fi

    ##### Mask-out data outside the brain ####
    #${FSLDIR}/bin/fslmaths "$T1wOutputDirectory"/data -mas "$T1wOutputDirectory"/nodif_brain_mask_temp "$T1wOutputDirectory"/data      ## edited from HCP script to add output checks
    if [ ! -f $T1wOutputDirectory/data_masked.nii.gz ] ; then
        echo "${i} Diffusion to Structural Registration : Masking-out data outside the brain"
        ${FSLDIR}/bin/fslmaths $T1wOutputDirectory/data -mas $T1wOutputDirectory/nodif_brain_mask_temp $T1wOutputDirectory/data_masked  
    else
        echo "${i} Diffusion to Structural Registration : Data Outside the Brain Masked-out"
    fi
    
    ##### Remove negative intensity values (caused by spline interpolation) from final data
    #${FSLDIR}/bin/fslmaths "$T1wOutputDirectory"/data -thr 0 "$T1wOutputDirectory"/data        ## edited from HCP script to add output checks
    if [ ! -f $T1wOutputDirectory/data_masked_thr.nii.gz ] ; then
        echo "${i} Diffusion to Structural Registration : Removing negative intensity values (caused by spline interpolation) from final data"
        ${FSLDIR}/bin/fslmaths $T1wOutputDirectory/data_masked -thr 0 $T1wOutputDirectory/data_masked_thr      
    else
        echo "${i} Diffusion to Structural Registration : Negative Intensity Values (caused by spline interpolation) Removed from Final Data"
    fi
#    ${FSLDIR}/bin/imrm $T1wOutputDirectory/nodif_brain_mask_temp

    if [ ! -f $T1wOutputDirectory/nodif_brain_mask_new.nii.gz ] ; then
        ${FSLDIR}/bin/fslmaths $T1wOutputDirectory/data_masked_thr -Tmean $T1wOutputDirectory/temp
        #${FSLDIR}/bin/immv $T1wOutputDirectory/nodif_brain_mask.nii.gz $T1wOutputDirectory/nodif_brain_mask_old.nii.gz     ## edited from HCP script to add output checks
        ${FSLDIR}/bin/imcp $T1wOutputDirectory/nodif_brain_mask.nii.gz $T1wOutputDirectory/nodif_brain_mask_old.nii.gz
        ${FSLDIR}/bin/fslmaths $T1wOutputDirectory/nodif_brain_mask_old.nii.gz -mas $T1wOutputDirectory/temp $T1wOutputDirectory/nodif_brain_mask_new
#        ${FSLDIR}/bin/imrm $T1wOutputDirectory/temp
   fi

done
















    
















