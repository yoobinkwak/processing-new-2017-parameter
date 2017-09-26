FSLDIR=/usr/share/fsl/5.0 
HCPDIR=Pipelines-master
HCPPIPEDIR_Global=${HCPPIPEDIR}/global/scripts

for i in $@
do
    preproc_func=${i}/Preprocess_Function
    if [ ! -d ${preproc_func} ]
    then
        mkdir ${preproc_func}
    fi

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
    if [ ! -f ${mc_matrix_dir1}/motion_matrix ]
    then
        bash ${HCPDIR}/fMRIVolume/scripts/edit_MotionCorrection.sh ${mc_dir1} ${i}/REST_MB4_LR_SBREF_0*/20* ${i}/REST_MB4_LR_SBREF_SBREF*/20* motion_corrected motion_regressors ${mc_matrix_dir1} motion_matrix MCFLIRT
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
    if [ ! -f ${mc_matrix_dir2}/motion_matrix ]
    then
        bash ${HCPDIR}/fMRIVolume/scripts/edit_MotionCorrection.sh ${mc_dir2} ${i}/REST_MB4_LR_SBREF_0*/20* ${i}/REST_MB4_LR_SBREF_SBREF*/20* motion_corrected motion_regressors ${mc_matrix_dir2} motion_matrix FLIRT
    else
        echo ${i} motion corrected with flirt
    fi
done


