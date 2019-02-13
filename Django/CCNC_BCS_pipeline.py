# python3
# kcho
#
# Dependencies
# ------------
# GPU with cuda support
# FSL, Freesurfer, Ants
# roiExtraction --> from thalamus_project kcho


from kcho import *
import glob
import datetime
from multiprocessing import Pool, current_process, Queue
import string
import shlex



def align_ACPC(inImg, outImg, template):
    ''' 
    ACPC align T1w or T2w image to 0.7mm MNI Template to create native volume space 
    '''

    imgDir = dirname(inImg)
    commands = '''robustfov \
        -i {ingImg} \
        -m {imgDir}/roi2full.mat \
        -r {imgDir}/robustroi.nii.gz;

    convert_xfm \
        -omat {imgDir}/full2roi.mat \
        -inverse {imgDir}/roi2full.mat;

    flirt \
        -interp spline \
        -in {imgDir}/robustroi.nii.gz \
        -ref {template} \
        -omat {imgDir}/roi2std.mat \
        -out {imgDir}/acpc_final.nii.gz \
        -searchrx -30 30 \
        -searchry -30 30 \
        -searchrz -30 30;

    convert_xfm \
        -omat {imgDir}/full2std.mat \
        -concat {imgDir}/roi2std.mat \
        {imgDir}/full2roi.mat;

    aff2rigid \
        {imgDir}/full2std.mat \
        {imgDir}/xfms_acpc.mat';

    applywarp \
        --rel \
        --interp=spline \
        -i {inImg} \
        -r {template} \
        --premat={imgDir}/xfms_acpc.mat \
        -o {ACPC}
    '''.format(inImg = inImg,
               template = template,
               imgDir = imgDir,
               ACPC = outImg)
    if not isfile(outImg):
        run(command)

def fnirt_based_masking(inImg, template, template2mm):
    betDir = join(dirname(inImg), 'FNIRTbased_bet')
    inImg_path_name_wo = join(dirname(inImg), 
                              basename(inImg).split('.')[0])
    inImg_brain = join(dirname(inImg),
                       basename(inImg).split('.')[0] + '_brain.nii.gz')
    inImg_brain_mask = join(dirname(inImg),
                       basename(inImg).split('.')[0] + '_brain_mask.nii.gz')
    template2mm_mask = join(dirname(template),
                            basename(template)+'_brain_mask.nii.gz')

    #template
    try:
        os.mkdir(betDir)
    except:
        pass

    command = '''flirt \
        -interp spline \
        -dof 12 \
        -in {inImg} \
        -ref {template} \
        -omat {betDir}/roughlin.mat \
        -out {inImg_path_name_wo}_to_MNI_roughlin.nii.gz \
        -nosearch;

    fnirt \
        --in={inImg}  \
        --ref={template2mm}  \
        --aff={betDir}/roughlin.mat  \
        --refmask={template2mm_mask}  \
        --fout={betDir}/str2standard.nii.gz  \
        --jout={betDir}/NonlinearRegJacobians.nii.gz  \
        --refout={betDir}/IntensityModulatedT1.nii.gz  \
        --iout={betDir}/{inImg_path_name_wo}_to_MNI_nonlin.nii.gz  \
        --logout={betDir}/NonlinearReg.txt  \
        --intout={betDir}/NonlinearIntensities.nii.gz  \
        --cout={betDir}/NonlinearReg.nii.gz  \
        --config=FNIRTConfig;

    {FSLDIR}/bin/applywarp \
        --rel \
        --interp=spline \
        --in={inImg} \
        --ref={template} \
        -w {betDir}/str2standard.nii.gz \
        --out={inImg_path_name_wo}_to_MNI_nonlin.nii.gz;

    {FSLDIR}/bin/invwarp \
         --ref={template2mm} \
         -w {betDir}/str2standard.nii.gz \
         -o {betDir}/standard2str.nii.gz;

    {FSLDIR}/bin/applywarp \
             --rel \
             --interp=nn \
             --in={template2mm_mask} \
             --ref={inImg} \
             -w {betDir}/standard2str.nii.gz \
             -o {inImg_brain_mask};

    {FSLDIR}/bin/fslmaths {inImg} \
             -mas {inImg_brain_mask} \
            {inImg_brain}'''.format(
                inImg = inImg,
                template = template,
                betDir = betDir,
                inImg_path_name_wo = inImg_path_name_wo,
                template2mm = template2mm,
                template2mm_mask = template2mm_mask,
                inImg_brain = inImg_brain,
                inImg_brain_mask = inImg_brain_mask)
    if not isfile(inImg_brain):
        run(command)

def reorient_to_mni(inImg, outImg):
    command = 'fslreorient2std {inImg} {outImg}'.format(
        inImg = inImg,
        outImg = outImg)
    if not isfile(outImg):
        run(command)

def is_list_single_item(file_list:list, required_file_name:str):
    '''
    If the list is composed of
    - one string : return string
    - more than one string : stop the script, ask the user to check
    - no item : stop the script, ask the user to check

    Used to check matching MRI directories

    '''
    file_number = len(file_list)

    if file_number > 1:
        sys.exit('There are more than one {} directory, please rename unused {} directory'.format(
            required_file_name))
    elif file_number == 0:
        #sys.exit('there are more than one {} directory, please rename unused {} directory'.format(
            #required_file_name))
        pass

    else:
        return file_list[0]


class CcncBcsSettings(CcncSettings):
    '''
    CCNC BCS settings
    '''
    def hcp_source(self):
        with open(self.hcp_setup_script, 'r') as f:
            lines = f.readlines()
        for line in lines:
            try:
                command = shlex.split(line)
                (key, _, value) = line.partition("=")
                key = key.split('export ')[1]
                #print('key:{}, value:{}'.format(key, value))
                os.environ[key] = value
                print(key, value)
            except:
                pass

    def hcp_run(self, command):
        new_command = 'source {}; {}'.format(self.hcp_setup_script, command)
        #new_command = 'bash {}'.format(command)
        run(new_command)

    def __init__(self, data_loc):
        super().__init__(data_loc)


        # HCP
        self.hcp_global_scripts = '/home/kangik/bin/HCPpipelines/global/scripts'
        self.hcp_setup_script = '/home/kangik/bin/HCPpipelines/Examples/Scripts/SetUpHCPPipeline.sh'
        #self.hcp_source()
        self.hcp_dir = '/home/kangik/bin/HCPpipelines'
        self.hcp_dti_dir = join(self.hcp_dir, 'DiffusionPreprocessing')

        self.subject = re.search('\w{3}_BCS\d{3}_\w{2,4}', data_loc).group(0)
        self.group = re.search('(\w{3})_BCS\d{3}_\w{2,4}', data_loc).group(1)

        # directory naming
        self.modality_retext_dict = {
            'T1':'^T1_\d+$',
            'T2':'^T2_\d+$',
            'REST':'REST_MB4_LR_SBREF_\d+$',
            'REST_LR_SBREF':'REST_MB4_LR_SBREF_SBREF_\d+$',
            'REST_LR':'REST_MB1_BLIP_LR_\d+$',
            'REST_RL':'REST_MB1_BLIP_RL_\d+$',
            'DTI_B1000':'^DTI_MB3_LR_B1000_\d+$',
            'DTI_B2000':'^DTI_MB3_LR_B2000_\d+$',
            'DTI_B3000':'^DTI_MB3_LR_B3000_\d+$',
            'DTI_BLIP_LR':'^DTI_BLIP_LR_\d+$',
            'DTI_BLIP_RL':'^DTI_BLIP_RL_\d+$'}

        self.modality_dcm_count_dict = {
            'T1':224,
            'T2':224,
            'REST':250,
            'REST_LR_SBREF':1,
            'REST_LR':3,
            'REST_RL':3,
            'DTI_BLIP_LR':7,
            'DTI_BLIP_RL':7,
            'DTI_B1000':21,
            'DTI_B2000':31,
            'DTI_B3000':65}

        self.timepoint = basename(data_loc)
        self.preproc_dir = join(data_loc, 'preprocessed')
        self.protocol_name = 'BCS'

        # find matching directories
        self.sub_dirs = [join(data_loc, x) for x in os.listdir(data_loc)]
        self.modality_name_dir_dict = self.check_and_get_modality_dirs()

    def check_and_get_modality_dirs(self):
        modality_dict = {}
        for modality_name, retext in self.modality_retext_dict.items():
            matching_dir_list = [x for x in self.sub_dirs if re.search(retext, basename(x), re.IGNORECASE)]
            modality_dir = is_list_single_item(matching_dir_list, modality_name)
            modality_dict[modality_name] = modality_dir

        #print(modality_dict)
        return modality_dict

    def check_dcm_numbers_in_modality_dirs(self):
        self.modality_check_dict = {}
        for modality_name, dcm_source in self.modality_name_dir_dict.items():
            dicom_target_number = self.modality_dcm_count_dict[modality_name]
            dicom_number = len([x for x in os.listdir(dcm_source) if re.search('dcm|ima', x, re.IGNORECASE)])

            if dicom_target_number == dicom_number:
                self.modality_check_dict[modality_name] = dicom_number
            else:
                self.modality_check_dict[modality_name] = dicom_number
                #print('{} - target number:{}, actual number:{}'.format(
                    #dcm_source, dicom_target_number, dicom_number))
        return self.modality_check_dict

    def convert_all_dcm_to_nifti(self):
        for modality_name, dcm_source in self.modality_name_dir_dict.items():
            target_dir = join(self.preproc_dir, modality_name)
            try:
                os.makedirs(target_dir)
            except:
                pass

            if modality_name+'.nii.gz' not in os.listdir(target_dir):
                convert_dcm_to_nifti(dcm_source, target_dir, modality_name) 


class CcncBcsDtiSettings(CcncBcsSettings, dtiSettings):
    '''
    CCNC BCS DTI pipeline
    22th Nov 2018
    '''
    def __init__(self, data_loc):
        super().__init__(data_loc)

        # dcm2nii conversion 
        self.dti_grappa_factor = 3
        self.dti_dir = join(self.preproc_dir, 'DTI')

    
        self.subject_dir = data_loc
        self.echo_spacing = 0.87
        self.dwi_data = join(self.dti_dir, 'data.nii.gz')
        self.blip_lr = join(self.preproc_dir, 'DTI_BLIP_LR', 'DTI_BLIP_LR.nii.gz')
        self.blip_rl = join(self.preproc_dir, 'DTI_BLIP_RL', 'DTI_BLIP_RL.nii.gz')
        self.blip_b0_merged = join(self.dti_dir, 'blip_b0_merged.nii.gz')

        # B0 extracted from BLIP images
        self.blip_lr_b0 = join(self.dti_dir, 'DTI_BLIP_LR_b0s.nii.gz')
        self.blip_rl_new = join(self.dti_dir, 'DTI_BLIP_RL_intensity_normcorr_to_LR.nii.gz')
        self.blip_rl_b0 = join(self.dti_dir, 'DTI_BLIP_RL_b0s.nii.gz')
        self.blip_rl_b0_new = join(self.dti_dir, 'DTI_BLIP_RL_b0s_intensity_normcorr_to_LR.nii.gz')

        self.blip_lr_bval = join(self.preproc_dir, 'DTI_BLIP_LR', 'DTI_BLIP_LR.bval')
        self.blip_rl_bval = join(self.preproc_dir, 'DTI_BLIP_RL', 'DTI_BLIP_RL.bval')


        self.bval = join(self.dti_dir, 'bvals')
        self.bvec = join(self.dti_dir, 'bvecs')

        # --------
        # topup
        # --------

        self.acqparams = join(self.dti_dir, 'acqparams.txt')
        self.topup_config_file=join(os.environ['FSLDIR'], 'etc/flirtsch/b02b0.cnf')
        self.topup_out = join(self.dti_dir, 'topup_out')
        self.dwi_topup = join(self.dti_dir, 'data_topup.nii.gz')


        # --------
        # Eddy
        # --------

        self.bet_f = 0.2 # bet f value initialise
        self.topup_mask = join(self.dti_dir, 'data_topup_brain_mask.nii.gz')
        self.nodif_brain = join(self.dti_dir, 'nodif_brain.nii.gz')
        self.nodif_brain_mask = join(self.dti_dir, 'nodif_brain_mask.nii.gz')

        self.eddy_index = join(self.dti_dir, 'eddy_index.txt')
        self.eddy_out = join(self.dti_dir, 'eddy_out')
        self.dwi_data_eddy_out = join(self.dti_dir, 'eddy_out.nii.gz')
        self.eddy_out_neg_rm = join(self.dti_dir, 'eddy_out_neg_rm.nii.gz')

        self.bvecs_eddy_out = join(self.dti_dir, 'eddy_out.eddy_rotated_bvecs')
        self.motion_rms = join(self.dti_dir, 'eddy_out.eddy_restricted_movement_rms')
        self.motion_rms_plot = join(self.dti_dir, 'eddy_out_motion.png')

        # --------
        # DTIFIT
        # --------
        self.FA = join(self.dti_dir, 'DTI_eddy_neg_rm_FA.nii.gz')
        self.MD = join(self.dti_dir, 'DTI_eddy_neg_rm_MD.nii.gz')
        self.RD = join(self.dti_dir, 'DTI_eddy_neg_rm_RD.nii.gz')

        # --------
        # Bedpostx
        # --------
        self.bedpostx_prep_dir = join(self.preproc_dir, 'DTI_preprocessed')
        self.bedpostx_dir = join(self.preproc_dir, 'DTI_preprocessed.bedpostX')

    def run_bet_topup(self):
        btr = fsl.BET(in_file = self.dwi_topup,
                      frac = self.bet_f,
                      out_file = self.nodif_brain,
                      mask = True)
        if not isfile(self.nodif_brain):
            btr.run()

    def merge_dtis(self):
        self.dwi_list = [join(self.preproc_dir, x, x+'.nii.gz') for x in ['DTI_B1000', 'DTI_B2000', 'DTI_B3000']]
        self.bval_list  = [join(self.preproc_dir, x, x+'.bval') for x in ['DTI_B1000', 'DTI_B2000', 'DTI_B3000']]
        self.bvec_list   = [join(self.preproc_dir, x, x+'.bvec') for x in ['DTI_B1000', 'DTI_B2000', 'DTI_B3000']]

        if not isfile(self.dwi_data):
            command = 'fslmerge -a {} {}'.format(self.dwi_data, ' '.join(self.dwi_list))
            run(command)
        np.savetxt(self.bval, np.concatenate([np.loadtxt(x) for x in self.bval_list]))
        np.savetxt(self.bvec, np.concatenate([np.loadtxt(x) for x in self.bvec_list], 1))
        
    def basic_dti_preprocessing(self):
        '''
        Normalization
        '''
        # Load BLIP b value table
        self.blip_lr_bval_array = np.loadtxt(self.blip_lr_bval)
        self.blip_rl_bval_array = np.loadtxt(self.blip_rl_bval)
        self.blip_lr_b0_num = len(self.blip_lr_bval_array[self.blip_lr_bval_array==0])
        self.blip_rl_b0_num = len(self.blip_rl_bval_array[self.blip_rl_bval_array==0])

        # Load BLIP images (LR, RL)
        blip_lr_img = nb.load(self.blip_lr)
        blip_rl_img = nb.load(self.blip_rl)
        blip_lr_data = blip_lr_img.get_data()
        blip_rl_data = blip_rl_img.get_data()

        # Load bvalue files
        blip_lr_bval = np.loadtxt(self.blip_lr_bval)
        blip_rl_bval = np.loadtxt(self.blip_rl_bval)

        # The order of 0 value in bvalue file
        # which represent the order (number) of b0 images in the blip lr/rl image
        # - BCS BLIP LR / RL images also include diffusion weighted i
        # - Could not add single B0 only volumes to the protocol back then.
        b0_locs_blip_lr = np.where(blip_lr_bval==0)[0]
        b0_locs_blip_rl = np.where(blip_rl_bval==0)[0]

        # Select b0 data arrays (volumes) from Blip lr data
        blip_lr_b0_array = blip_lr_data[:,:,:,b0_locs_blip_lr]
        blip_rl_b0_array = blip_rl_data[:,:,:,b0_locs_blip_rl]

        # Average of the 4D b0 volumes
        blip_lr_b0_avg = np.mean(blip_lr_b0_array)
        blip_rl_b0_avg = np.mean(blip_rl_b0_array)

        # Normalize blip_rl using the average from blip_rl
        # HCP protocol
        # ${FSLDIR}/bin/fslmaths ${basename} -mul ${rescale} -div ${scaleS} ${basename}_new
        rescale = blip_lr_b0_avg
        scale = blip_rl_b0_avg
        blip_rl_data_new = blip_rl_data * rescale / scale
        blip_rl_b0_array_new = blip_rl_data_new[:,:,:,b0_locs_blip_rl]

        # merge blip LR + RL b0 images
        blip_b0_merged_data = np.concatenate([blip_lr_b0_array, blip_rl_b0_array_new], 3)

        # Save arrays into nifti files
        # - blip_lr_b0
        # - blip_rl_new - normalized blip_rl
        # - blip_rl_b0_new - normalized blip_rl_b0
        for img_loc, array in zip([self.blip_lr_b0,
                                   self.blip_rl_new,
                                   self.blip_rl_b0_new,
                                   self.blip_b0_merged],
                                  [blip_lr_b0_array,
                                   blip_rl_data_new,
                                   blip_rl_b0_array_new,
                                   blip_b0_merged_data]):
            if not isfile(img_loc):
                nb.Nifti1Image(array, 
                               affine=blip_lr_img.affine).to_filename(img_loc)


    def write_acqparams(self):
        '''
        A >> P : 0 -1 0
        A << P : 0 1 0
        4th number : 0.69 ms * 112 * 0.001
        https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind1306&L=fsl&D=0&P=49368
        > TR = 9300ms , TE=94ms, Echo spacing = 0.69ms, 96x96 matrix and 65 slices,
        > Phase partial Fourier 6/8 and finally bandwidth 1628Hz/Px
        the relevant time in your case is (96-1)*0.00069 = 0.0656 seconds.
        For SCS project
        - 112x112 matrix
        - 0.8 echo spacing
        (112-1) * 0.0008 = .0888
        '''
        #num = (matrix_size-1) * (1/1000 * echo_spacing)
        #acqparamLoc = join(outDir, 'acqparams.txt')

        #if not isfile(acqparamLoc):
        if not isfile(self.acqparams):
            # Writing acqparams.txt
            lr_array = np.tile([1, 0, 0, self.ro_time], self.blip_lr_b0_num).reshape(self.blip_lr_b0_num, 4)
            rl_array = np.tile([-1, 0, 0, self.ro_time], self.blip_rl_b0_num).reshape(self.blip_rl_b0_num, 4)
            concat_array = np.concatenate([lr_array, rl_array])
            np.savetxt(self.acqparams, concat_array, 
                       fmt=['%d', '%d', '%d', '%0.3f'])


    def run_topup(self):
        command = 'topup \
                --imain={blip_b0_merged} \
                --datain={acqparams} \
                --config={topup_config_file} \
                --out={topup_out} \
                -v'.format(blip_b0_merged=self.blip_b0_merged,
                        acqparams = self.acqparams,
                        topup_config_file = self.topup_config_file,
                        topup_out = self.topup_out)
        if not isfile(self.topup_out + '_fieldcoef.nii.gz'):
            #print(command)
            run(command)

    def apply_topup(self):
        command = 'applytopup \
                --imain={b0_lr},{b0_rl} \
                --topup={topup_out} \
                --datain={acqparams} \
                --inindex={b0_lr_num},{b0_rl_num} \
                --out={dwi_data_topup}'.format(
                        b0_lr=self.blip_lr_b0,
                        b0_rl=self.blip_rl_b0_new,
                        acqparams = self.acqparams,
                        b0_lr_num = 1,
                        b0_rl_num = 7,
                        topup_out = self.topup_out,
                        dwi_data_topup = self.dwi_topup)

        if not isfile(self.dwi_topup):
            run(command)

    def eddy(self):
        # index
        bvals_mat = np.loadtxt(self.bval)
        vol_num = bvals_mat.shape[0]

        # create an index file
        # 70 --> number of volumes SCS project
        index = np.tile(1, vol_num)
        np.savetxt(self.eddy_index, index, fmt='%d', newline=' ')

        #eddy
        #command = 'eddy_cuda8.0 \
        command = 'eddy \
                --imain={dwi} \
                --mask={mask} \
                --acqp={acq} \
                --index={index_loc} \
                --bvecs={bvecs} \
                --bvals={bvals} \
                --fwhm=0 \
                --flm=quadratic \
                --topup={topup_prefix} \
                --out={eddy_out}'.format(
                        dwi = self.dwi_data,
                        mask = self.nodif_brain_mask,
                        acq = self.acqparams,
                        index_loc = self.eddy_index,
                        bvecs = self.bvec,
                        bvals = self.bval,
                        topup_prefix = self.topup_out,
                        eddy_out = self.eddy_out)

        
        #gpu_command = 'CUDA_VISIBLE_DEVICES={} {}'.format(gpu_num, command)

        if not isfile(self.eddy_out+'.nii.gz'):
            run(command)
            #run_gpu_command(command, gpu_num)
            #print(re.sub('\s+', ' ', command))

    def copy_eddy_out_to_data(self):
        shutil.copy(self.eddy_out_neg_rm, self.dti_dir, 'data.nii.gz')

    def eddy_remove_negative_intensities(self):
        command = 'fslmaths {eddy_out} \
                -thr 0 \
                {eddy_out_neg_removed}'.format(
                    eddy_out = self.eddy_out,
                    eddy_out_neg_removed = self.eddy_out_neg_rm)
        if not isfile(self.eddy_out_neg_rm):
            run(command)

    def dtifit_eddy_out(self):
        command = 'dtifit -k {data} -o {outname} -m {mask} -r {bvecs} -b {bvals}'.format(
            data = self.eddy_out_neg_rm,
            outname = join(self.dti_dir, 'DTI_eddy_neg_rm'),
            mask = self.nodif_brain_mask,
            bvecs = self.bvecs_eddy_out,
            bvals = self.bval)

        if not isfile(self.FA):
            run(command)

    def bedpostx(self):
        # copy uddy processed files
        if not isfile(join(self.bedpostx_prep_dir, 'nodif_brain_mask.nii.gz')):
            try:
                os.mkdir(self.bedpostx_prep_dir)
            except:
                pass
            for initial_name, new_name in zip([self.bval, 
                                               self.bvec, 
                                               self.eddy_out_neg_rm, 
                                               self.nodif_brain_mask],
                                              ['bvals', 
                                               'bvecs', 
                                               'data.nii.gz',
                                               'nodif_brain_mask.nii.gz']):
                shutil.copy(initial_name,
                            join(self.bedpostx_prep_dir, new_name))

        if not isdir(self.bedpostx_dir):
            command = 'bedpostx {dtiDir} -n 3 -model 2'.format(dtiDir = self.bedpostx_prep_dir)
            #gpu_command = 'CUDA_VISIBLE_DEVICES={} {}'.format(gpu_num, command)
            #print(gpu_command)
            run(command)

    def draw_motion_graph(self):
        '''
        motion_array : (n,2) matrix
        absolute and relative motion parameters from eddy
        '''
        if not isfile(self.motion_rms_plot):
            motion_array = np.loadtxt(self.motion_rms)

            fig, axes = plt.subplots(nrows=2, dpi=150)

            absolute_ax = axes[0]
            relative_ax = axes[1]

            for (num, ax), name, color in zip(enumerate(np.ravel(axes)),
                                              ['Absolute', 'Relative'],
                                              ['b', 'r']):
                ax.set_title(name + ' motion')
                ax.plot(motion_array[:,num], 
                        color=color, 
                        label = name + ' motion')
                ax.axhline(2, linestyle='--')
                ax.set_ylim(0, 3)
                ax.set_xlim(-5, 120)
                ax.set_ylabel('Displacement in mm')
                ax.set_xlabel('Diffusion data timepoint')

            fig.suptitle('Motion in diffusion : {} {}'.format(self.subject, 
                                                              self.timepoint))
            fig.tight_layout()
            plt.subplots_adjust(top=0.88)
            fig.savefig(self.motion_rms_plot)
        else:
            pass
        

class CcncBcsStructSettings(CcncBcsSettings):
    def __init__(self, data_loc):
        super().__init__(data_loc)
        self.fs_dir = join(self.preproc_dir, 'FREESURFER')
        
        self.t1_dir = join(self.preproc_dir, 'T1')
        self.t2_dir = join(self.preproc_dir, 'T2')
        self.t1 = join(self.preproc_dir, 'T1/T1.nii.gz')
        self.t2 = join(self.preproc_dir, 'T2/T2.nii.gz')

        self.struct_preproc_dir = join(self.preproc_dir, 'Struct_preproc')
        self.t1_acpc_dc = join(self.struct_preproc_dir, 'T1w/T1w_acpc_dc_restore.nii.gz')
        self.t1_acpc_dc_brain = join(self.struct_preproc_dir, 'T1w/T1w_acpc_dc_restore_brain.nii.gz')
        self.brainmask_fs = join(self.struct_preproc_dir, 'T1w/T1w_acpc_brain_mask.nii.gz') 
        self.t1_acpc_dc_biasfield = join(self.struct_preproc_dir, 'T1w/BiasField_acpc_dc.nii.gz') 
        self.t2_acpc_dc = join(self.struct_preproc_dir, 'T1w/T2w_acpc_dc_restore.nii.gz')



        self.fs_dir = join(self.preproc_dir, 'FREESURFER')

    def run_hcp_freesurfer(self):
        command = '{script} \
                    --subject={subject} \
                    --subjectDIR={SUBJECTS_DIR} \
                    --t1={t1_acpc_dc} \
                    --t1brain={t1_acpc_dc_brain} \
                    --t2={t2_acpc_dc}'.format(
                        script = self.hcp_fs_script,
                        subject = 'FREESURFER',
                        SUBJECTS_DIR = self.preproc_dir,
                        t1_acpc_dc = self.t1_acpc_dc,
                        t1_acpc_dc_brain = self.t1_acpc_dc_brain,
                        t2_acpc_dc = self.t2_acpc_dc)
        print(command)
        self.hcp_run(command)


    def run_hcp_prefreesurfer(self):
        command = '{script} \
                --path={PATH} \
                --subject={SUBJ} \
                --t1={T1} \
                --t2={T2}'.format(
                    script = self.hcp_pre_fs_script,
                    PATH = self.preproc_dir,
                    SUBJ = 'Struct_preproc',
                    T1 = self.t1,
                    T2 = self.t2)

        if not isfile(self.t1_acpc_dc_brain):
        #if not isfile(self.t2_acpc_dc):
            self.hcp_run(command)

class CcncBcsStructPreproc(CcncBcsStructSettings):
    def __init__(self, data_loc):
        super().__init__(data_loc)
        try:
            self.run_hcp_prefreesurfer()
            self.run_hcp_freesurfer()
        except:
            pass


class CcncBcsDtiPreproc(CcncBcsDtiSettings):
    def __init__(self, data_loc):
        super().__init__(data_loc)

        #try:
        # dcm2nii conversion 
        pheadline('Convert all dicoms into niftis')
        self.check_dcm_numbers_in_modality_dirs()
        self.convert_all_dcm_to_nifti()

        # make DTI directory
        pheadline('Make DTI preproc directory')
        try:
            os.mkdir(self.dti_dir)
        except:
            pass
    
        # merge DTIs
        pheadline('Merge B1000, B2000, B3000 DTIs')
        self.merge_dtis()

        # normalize blip rl with blip lr
        # save extracted and merged b0 volumes
        pheadline('Basic preprocessing')
        self.basic_dti_preprocessing()

        # Compute Total_readout in secs with up to 6 decimal places 
        # self.dimP=`${FSLDIR}/bin/fslval ${rawdir}/RL_BLIP_input.nii.gz dim1`
        # nPEsteps=$(($dimP - 1))
        # ro_time=`echo "${echo_spacing} * (${nPEsteps}/${grappa_factor})" | bc -l`
        # ro_time=`echo "scale=6; ${ro_time} / 1000" | bc -l`
        self.grappa_factor = 3
        self.dimP = nb.load(self.blip_lr).get_data().shape[0]
        self.nPEsteps = self.dimP - 1
        self.ro_time = (self.echo_spacing * (self.nPEsteps / self.grappa_factor)) / 1000

        # Topup
        # Takes about 1hr in CCNC GPU server
        pheadline('Topup')
        self.write_acqparams()
        self.run_topup()
        self.apply_topup()

        ## Eddy
        # Takes about 25 mins in CCNC GPU server
        pheadline('Eddy')
        self.run_bet_topup()
        self.eddy()
        self.eddy_remove_negative_intensities()
        #self.draw_motion_graph()

        ## dtifit
        self.dtifit_eddy_out()

        # bedpostx
        self.bedpostx()
        #try:
            ##pheadline('Error in subject '+ data_loc)
            ### edit here later
            ##self.eddy_remove_negative_intensities()
            ###self.draw_motion_graph()
            ##self.dtifit_eddy_out()
            ### Eddy
            ## Takes about 25 mins in CCNC GPU server
            #pheadline('Eddy')
            #self.run_bet_topup()
            #self.eddy()
            #self.eddy_remove_negative_intensities()
            #self.copy_eddy_out_to_data()
            ##self.draw_motion_graph()

            ### dtifit
            #self.dtifit_eddy_out()

            ## bedpostx
            ##self.bedpostx()
        #except:
            #pheadline('Error in subject '+ data_loc)
            ## edit here later
            #self.eddy_remove_negative_intensities()
            #self.copy_eddy_out_to_data()
            ##self.draw_motion_graph()
            #self.dtifit_eddy_out()


class CcncBcsDatabase(CcncBcsDtiSettings):
    '''
    Update database
    '''
    def __init__(self, data_root):
        self.data_root = data_root

        self.modality_dcm_count_dict = {
            'T1':224,
            'T2':224,
            'REST':250,
            'REST_LR_SBREF':1,
            'REST_LR':3,
            'REST_RL':3,
            'DTI_BLIP_LR':7,
            'DTI_BLIP_RL':7,
            'DTI_B1000':21,
            'DTI_B2000':31,
            'DTI_B3000':65}
        
        # Make data frame using glob.glob : creates self.df
        self.get_basic_df_from_data_points()

        # Add demographics from the MRI team DB : updates self.df
        self.db_loc_orig = join(data_root, 'database/database_bcs.xlsx')
        self.update_from_prev_db()

        # Count dicoms : updates self.df
        self.make_dcm_count_df()

        self.df_loc = join(data_root, 'prac.xlsx')
        self.writer = pd.ExcelWriter(self.df_loc, engine='xlsxwriter', 
                                     datetime_format='yyyy-mm-dd',
                                     date_format='yyyy-mm-dd')

        # Update DTI preprocessing : updates self.df
        self.make_preprocessing_df()

        # DTI head motion extraction
        #self.get_head_motion_DTI()

        # All data df write
        self.df.to_excel(self.writer, sheet_name='Raw data check')

        # Format excelfile
        self.format_excel()

        # DTI preproc sheet write
        #self.dti_df.to_excel(self.writer, sheet_name='DTI preprocessing')
        #worksheet = self.writer.sheets['DTI preprocessing']
        self.check_counts()
        self.writer.save()

        #self.run_missing_topup()

        # Preprocess T1
        #for i, row in self.df.groupby('T1')
        #CcncBcsDtiPreproc(self.df.loc[0, 'Data point'])
        print(self.df.loc[0, 'Data points'])
        #CcncBcsStructPreproc(self.df.loc[0, 'Data points'])

        #CcncBcsStructPreproc(self.df.loc[0, 'Data points'])
        CcncBcsDtiStructPreproc('/Volumes/CCNC_4T/BCS_MRI/CHR/CHR_BCS001_BBS/baseline')
        #CcncBcsDtiStructPreproc(self.df.loc[0, 'Data points'])

        #with Pool(procEsses=12) as pool:
           #results = pool.map(CcncBcsStructPreproc, self.df['Data points'])

    def get_head_motion_DTI(self):
        self.dti_motion_df = self.df.groupby('Eddy').get_group(True)
        self.dti_motion_df['Eddy motion file'] = self.dti_motion_df['Data points'] + '/preprocessed/DTI/eddy_out.eddy_restricted_movement_rms'
        self.dti_motion_df['Eddy motion array'] = self.dti_motion_df['Eddy motion file'].apply(lambda x: np.loadtxt(x))
        self.dti_motion_df['Eddy motion absolute average'] = self.dti_motion_df['Eddy motion array'].apply(lambda x: np.mean(x[:,0]))
        self.dti_motion_df['Eddy motion relative average'] = self.dti_motion_df['Eddy motion array'].apply(lambda x: np.mean(x[:,1]))
        self.dti_motion_df['Eddy motion absolute std'] = self.dti_motion_df['Eddy motion array'].apply(lambda x: np.std(x[:,0]))
        self.dti_motion_df['Eddy motion relative std'] = self.dti_motion_df['Eddy motion array'].apply(lambda x: np.std(x[:,1]))

        # Absolute any point over 2mm
        fig, ax = plt.subplots(nrows=1, figsize=(10,5))
        for index, row in self.dti_motion_df.iterrows():
            if np.any(row['Eddy motion array'][:,0]>2): 
                ax.plot(row['Eddy motion array'][:,0], alpha=0.6, label=row['Subject'][:10])
                print(row['Subject'][:10])

        ax.axhline(2, linestyle='--')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        fig.legend(loc = 'center left', bbox_to_anchor=(1, 0.5))
        #fig.legend()
        fig.tight_layout()
        fig.savefig('prac2.png', bbox_inches='tight')

        # Relative any point over 2mm
        fig, ax = plt.subplots(nrows=1, figsize=(10,5))
        for index, row in self.dti_motion_df.iterrows():
            if np.any(row['Eddy motion array'][:,1]>2): 
                ax.plot(row['Eddy motion array'][:,1], alpha=0.6, label=row['Subject'][:10])
                print(row['Subject'][:10])

        ax.axhline(2, linestyle='--')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        fig.legend(loc = 'center left', bbox_to_anchor=(1, 0.5))
        #fig.legend()
        fig.tight_layout()
        fig.savefig('prac3.png', bbox_inches='tight')
        
        self.df['DTI absolute motion average'] = self.dti_motion_df['Eddy motion absolute average']
        self.df['DTI absolute motion std'] = self.dti_motion_df['Eddy motion absolute std']
        self.df['DTI relative motion average'] = self.dti_motion_df['Eddy motion relative average']
        self.df['DTI relative motion std'] = self.dti_motion_df['Eddy motion relative std']

        melt_df = pd.melt(self.dti_motion_df, id_vars=['Subject', 'Timepoint', 'Group'], 
                          value_vars=['Eddy motion absolute average', 'Eddy motion relative average'], 
                          var_name='Motion type', value_name='Motion')
        g = sns.catplot(x='Subject', y='Motion', hue='Group', col='Motion type',
                       data=melt_df)#.savefig(join(self.data_root, "prac.png"))

        #for ax in g.ax:
        for ax, motion_type in zip(np.ravel(g.axes), ['Eddy motion absolute average', 'Eddy motion relative average']):
            ax.axhline(2, linestyle='--')
            ax.get_xaxis().set_ticks([])

            ax.set_title(motion_type)
            # subjects with high motion
            motion_threshold = 2
            tmp_df = melt_df.groupby('Motion type').get_group(motion_type)
            x_indicies = tmp_df[tmp_df['Motion'] >= 2].index
            y_indicies = tmp_df[tmp_df['Motion'] >= 2]['Motion']
            subject_names = tmp_df[tmp_df['Motion'] >= 2]['Subject'].str[:10]

            for x, y, subject_name in zip(x_indicies, y_indicies, subject_names):
                ax.text(x, y+0.1, subject_name, horizontalalignment='center')

        g.savefig(join(self.data_root, "prac.png"))

        # Group motion
        df = pd.DataFrame()
        for group in self.dti_motion_df['Group'].unique():
            group_df = self.dti_motion_df.groupby('Group').get_group(group)
            group_average_motion_ts = np.mean(np.stack(group_df['Eddy motion array'], axis = -1), 2)
            df_tmp = pd.DataFrame(group_average_motion_ts, columns=['Group absolute motion', 
                                                                    'Group relative motion'])
            df_tmp['Group'] = group
            df = pd.concat([df, df_tmp])

        self.dti_motion_df = df

        # Plot group motion
        fig, axes = plt.subplots(nrows=2, figsize=(10,7))
        for group in self.dti_motion_df['Group'].unique():
            df = self.dti_motion_df.groupby('Group').get_group(group)
            axes[0].plot(df['Group absolute motion'], label=group)
            axes[0].set_title('Group absolute motion')

            axes[1].plot(df['Group relative motion'], label='_nolegend_')
            axes[1].set_title('Group relative motion')
            axes[1].set_xlabel('Diffusion data timepoint')

        for ax in np.ravel(axes):
            ax.axhline(2, linestyle='--')
            ax.set_ylim(0, 3)
            ax.set_xlim(-5, 120)
            ax.set_ylabel('Displacement in mm')
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        fig.suptitle('DTI motion for each group', x=0.53)
        fig.legend(loc = 'center left', bbox_to_anchor=(1, 0.5))
        fig.tight_layout()
        plt.subplots_adjust(top=0.90)
        fig.savefig('./prac.png', bbox_inches='tight')

    def update_from_prev_db(self):
        db_df = pd.read_excel(self.db_loc_orig, 'Sheet1')
        to_select_from_db = ['koreanName', 'sex', 'age', 'DOB', 'scanDate', 'timeline',
                             'studyname', 'patientNumber', 'folderName', 'dx', 'note']
        db_df = db_df[to_select_from_db]

        #to_select_from_data = ['Subject', 'Timepoint', 'Group',
                               #'T1', 'T2',
                               #'REST', 
                               #'REST_LR_SBREF', 'REST_LR', 'REST_RL',
                               #'DTI_BLIP_LR', 'DTI_BLIP_RL', 
                               #'DTI_B1000','DTI_B2000','DTI_B3000']
        to_select_from_data = ['Subject', 'Timepoint', 'Group',
                               'T1', 'T2',
                               'REST', 
                               'REST_LR_SBREF', 'REST_LR', 'REST_RL',
                               'DTI_BLIP_LR', 'DTI_BLIP_RL', 
                               'DTI_B1000','DTI_B2000','DTI_B3000']
        #tmp_df = self.real_time_db#[to_select_from_data]
        tmp_df = pd.merge(db_df, self.df,
                          left_on=['folderName', 'timeline'],
                          right_on=['Subject', 'Timepoint'],
                          how='outer')#.drop('timeline', axis=1).drop('folderName', axis=1)
        self.df = tmp_df

    def make_dcm_count_df(self):
        '''
        dicom count sheet
        '''

        # Get modality and dicom count dictionary 
        self.df['tmp_dict'] = self.df['Data points'].apply(
            lambda x: CcncBcsSettings(x).check_dcm_numbers_in_modality_dirs() if type(x)==str else 0)

        for modality in self.modality_dcm_count_dict.keys():
            self.df[modality] = self.df['tmp_dict'].str[modality]
            target_dicom_number = self.modality_dcm_count_dict[modality]
            #df[modality+'_check'] = df[modality]==target_dicom_number

        columns = ['koreanName', 'Subject', 'Timepoint', 'Group', 
                   'sex', 'age', 'DOB', 'scanDate', 'studyname', 'patientNumber', 'dx',
                   'T1', 'T2', 'REST', 'REST_LR_SBREF',
                   'REST_LR', 'REST_RL', 'DTI_BLIP_LR', 'DTI_BLIP_RL', 'DTI_B1000',
                   'DTI_B2000', 'DTI_B3000',
                   'note', 'Data points']

        self.df = self.df[columns]

    def check_counts(self):
        # Group x Timepoint count
        save_count = lambda gb, name: gb['Subject'].agg('count').to_frame('count').to_excel(self.writer, sheet_name=name)
        save_count(self.df.groupby(['Group', 'Timepoint']), 'Subject Count')

        modality_group_dict = {'T1':['T1'], 
                               'T2':['T2'],
                               'DTI':['DTI_BLIP_LR', 'DTI_BLIP_RL', 'DTI_B1000', 'DTI_B2000', 'DTI_B3000'],
                               'REST':['REST','REST_LR_SBREF', 'REST_LR', 'REST_RL']}

        tmp_df = self.df.copy()
        # for T1, T2, DTI and REST modalities
        for modality_group, modality_groups in modality_group_dict.items():
            # For each modalities for each modality group
            for modality in modality_groups:
                target_dicom_number = self.modality_dcm_count_dict[modality]
                # check whether the number of dicoms match the target number
                tmp_df[modality+'_check'] = (tmp_df[modality] == target_dicom_number).map({True:'with', 
                                                                           False:'without'})

            # Make a new column checking all modality groups
            tmp_df[modality_group] = np.all([tmp_df[modality + '_check'] for modality in modality_groups], axis=0)
            save_count(tmp_df.groupby(['Group', 'Timepoint', modality_group]), modality_group)

    def make_preprocessing_df(self):
        #df = self.df.copy()

        # Data preprocessing check
        self.df['Preprocessing started'] = self.df['Data points'].apply(lambda x: isdir(join(x, 'preprocessed')) if type(x)==str else False)

        # DTI
        column_file_dict = {'Basic preprocessing':'data.nii.gz',
                            'Topup': 'topup_out_fieldcoef.nii.gz',
                            'Apply Topup' : 'data_topup.nii.gz',
                            'Bet for Eddy' : 'nodif_brain_mask.nii.gz',
                            'Eddy' : 'eddy_out.nii.gz',
                            'Eddy rm neg' : 'eddy_out_neg_rm.nii.gz',
                            'DTIFIT' : 'DTI_eddy_neg_rm_FA.nii.gz',
                            'Bedpostx' : ''}

        for col, file_name in column_file_dict.items():
            self.df[col] = self.df['Data points'].apply(lambda x: isfile(join(x, 'preprocessed/DTI', file_name) if type(x)==str else False))

        # T1
        #column_file_dict = {'
        self.df['Bedpostx'] = self.df['Data points'].apply(lambda x: isfile(join(x, 'preprocessed', 'DTI_preprocessed.bedpostX', 'mean_fsumsamples.nii.gz')) if type(x)==str else False)
        self.df['DTI preprocessing'] = np.all([self.df[stage] for stage in column_file_dict.keys()], axis=0)

        self.dti_df = self.df[['Subject', 'Timepoint', 'Group'] + \
                          list(column_file_dict.keys()) + \
                          ['DTI preprocessing', 'Data points']]

    def format_excel(self):

        # Excel file settings : Light red fill with dark red text.
        worksheet = self.writer.sheets['Raw data check']
        workbook = self.writer.book
        self.format1 = workbook.add_format({'bg_color':   '#FFC7CE',
                                            'font_color': '#9C0006'})
        ## Get column order for T1
        #T1_column_order = self.df.columns.index('T1')
        #number_to_alphabet_dict=dict(zip(string.ascii_uppercase,[ord(c)%32 for c in string.ascii_uppercase]))
        #alphabet_for_t1_column = number_to_alphabet_dict[dicom_count_column_order+1]

        # dicom numbers
        for num, modality in enumerate(self.modality_dcm_count_dict.keys()):
            alphabet = 'MNOPQRSTUVW'[num]
            target_dicom_number = self.modality_dcm_count_dict[modality]
            worksheet.conditional_format('{col}2:{col}{length}'.format(col = alphabet, length = len(self.df)+1),
                                         {'type': 'cell',
                                          'criteria': '!=',
                                          'value': target_dicom_number,
                                          'format': self.format1})

        # preprocessing steps
        worksheet.conditional_format('Z2:AI{length}'.format(length = len(self.df)+1),
                                     {'type': 'cell',
                                      'criteria': '!=',
                                      'value': 'True',
                                      'format': self.format1})
        # DTI motion format
        worksheet.conditional_format('AJ2:AM{length}'.format(length = len(self.df)+1),
                                     {'type': 'cell',
                                      'criteria': '>',
                                      'value': 2,
                                      'format': self.format1})
        # Column width
        worksheet.set_column('A:AI', 18)


    def get_basic_df_from_data_points(self):
        # Catch data points
        self.data_point_list = [x for x in glob.glob(self.data_root + '/*/*/*') \
                                if re.search('\/[A-Z]{3}/[A-Z]{3}_BCS\d{3}_\w{2,4}', x)]

        # Folder directory db
        self.df = pd.DataFrame(self.data_point_list, columns=['Data points'])
        self.df['Subject'] = self.df['Data points'].str.extract('\/([A-Z]{3}_BCS\d{3}_\w{2,4})\/')
        self.df['Timepoint'] = self.df['Data points'].apply(basename)
        self.df['Group'] = self.df['Subject'].str.extract('([A-Z]{3})')

    def run_missing_topup(self):
        #data_locs = self.dti_df.groupby('Topup').get_group(False)['Data points']
        #print(data_locs)
        #data_locs = self.dti_df['Data points'].unique()
        #with Pool(processes=12) as pool:
            #pheadline('Running DTI preprocessing upto Bet in {} subjects'.format(len(data_locs)))
            #results = pool.map(CcncBcsDtiPreproc, data_locs)
        
        # tmp remove below later
        # eddy run
            #tmp = CcncBcsDtiSettings(data_loc)
            ## Eddy
            #tmp.run_bet_topup()
            #tmp.eddy()
            #tmp.eddy_remove_negative_intensities()
            #tmp.draw_motion_graph()

            # dtifit
            #tmp.dtifit_eddy_out()

            #bedpostx
            #tmp.bedpostx()

        NUM_GPUS = 7
        PROC_PER_GPU = 1

        #global queue

        data_locs_eddy_to_run = self.dti_df.groupby(['Bedpostx']).get_group(False)['Data points']
        print(data_locs_eddy_to_run)

        # initialize the queue with the GPU ids
        for gpu_ids in range(NUM_GPUS):
            for _ in range(PROC_PER_GPU):
                queue.put(gpu_ids)

        pool = Pool(processes=PROC_PER_GPU * NUM_GPUS)
        for _ in pool.imap_unordered(foo, data_locs_eddy_to_run):
            pass
        pool.close()
        pool.join()


class CcncBcsDtiRegSettings(CcncBcsDtiSettings, CcncBcsStructSettings):
    def __init__(self, data_loc):
        super().__init__(data_loc)
        self.dti_struct_reg_dir = join(self.preproc_dir, 'registration')
        self.dti_in_struct_dir = join(self.preproc_dir, 'DTI_in_struct_space')

        self.dti_struct_reg_tmp = join(self.preproc_dir, 
                                       'registration',
                                       'reg_tmp')

        self.dti_struct_qa_tmp = join(self.preproc_dir, 
                                       'registration',
                                       'qa_tmp')


    def dti_post_eddy_processing(self):
        try:
            os.mkdir(self.dti_struct_reg_dir)
        except:
            pass

        try:
            os.mkdir(self.dti_in_struct_dir)
        except:
            pass

        command = '{hcp_dti_dir}/scripts/DiffusionToStructural.sh \
                --t1folder={t1_dir} \
                --subject=FREESURFER \
                --workingdir={reg_dir} \
                --datadiffdir={dti_dir} \
                --t1={t1} \
                --t1restore={t1_restore} \
                --t1restorebrain={t1_restore_brain} \
                --biasfield={biasfield} \
                --brainmask={fs_brain_mask} \
                --diffresol=2.3 \
                --datadiffT1wdir={dti_in_struct_dir} \
                --regoutput={tmp_regout_file} \
                --QAimage={tmp_qa_file}'.format(
                    hcp_dti_dir = self.hcp_dti_dir,
                    t1_dir = self.preproc_dir,
                    reg_dir = self.dti_struct_reg_dir,
                    dti_dir = self.dti_dir,
                    t1 = self.t1.split('.nii.gz')[0],
                    t1_restore = self.t1_acpc_dc.split('.nii.gz')[0],
                    t1_restore_brain = self.t1_acpc_dc_brain.split('.nii.gz')[0],
                    biasfield = self.t1_acpc_dc_biasfield,
                    fs_brain_mask = self.brainmask_fs,
                    dti_in_struct_dir = self.dti_in_struct_dir,
                    tmp_regout_file = self.dti_struct_reg_tmp,
                    tmp_qa_file = self.dti_struct_qa_tmp)
        print(command)

class CcncBcsDtiStructPreproc(CcncBcsDtiRegSettings):
    def __init__(self, data_loc):
        super().__init__(data_loc)

        # run post freesurfer * post eddy processing
        # HCP
        self.dti_post_eddy_processing()

queue = Queue()
def foo(data_loc):
    gpu_id = queue.get()
    print(gpu_id)
    try:
        # run processing on GPU <gpu_id>
        ident = current_process().ident
        print('{}: starting process on GPU {}'.format(ident, gpu_id))
        tmp = CcncBcsDtiSettings(data_loc)
        #tmp.run_bet_topup()
        #tmp.eddy(gpu_id)
        ha = tmp.bedpostx(gpu_id)
        print(ha)
        # ... process filename
        print('{}: finished'.format(ident))
    except:
        pass
    finally:
        queue.put(gpu_id)


if __name__ == "__main__":
    data_root = '/Volumes/CCNC_4T/BCS_MRI'
    tmp = CcncBcsDatabase(data_root)
    #now = datetime.datetime.now()
    #self.datetime = '{}_{}_{}'.format(now.year, now.month, now.day)

