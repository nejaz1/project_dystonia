
mri_convert d01_anatomical.nii d01_anatomical.mgh

setenv SUBJECTS_DIR ~/Projects/FingerPattern_dystonia/surfaceFreesurfer/
cd ~/Projects/FingerPattern_dystonia/surfaceFreesurfer/
recon-all -subject d01 -i d01_anatomical.mgh -all -cw256


mri_convert s01_anatomical.nii s01_anatomical.mgh

setenv SUBJECTS_DIR ~/Projects/FingerPattern_dystonia/surfaceFreesurfer/
cd ~/Projects/FingerPattern_dystonia/surfaceFreesurfer/
recon-all -subject s01 -i s01_anatomical.mgh -all -cw256



mri_convert d02_anatomical.nii d02_anatomical.mgh

setenv SUBJECTS_DIR ~/Projects/FingerPattern_dystonia/surfaceFreesurfer/
cd ~/Projects/FingerPattern_dystonia/surfaceFreesurfer/
recon-all -subject d02 -i d02_anatomical.mgh -all -cw256


mri_convert s02_anatomical.nii s02_anatomical.mgh

setenv SUBJECTS_DIR ~/Projects/FingerPattern_dystonia/surfaceFreesurfer/
cd ~/Projects/FingerPattern_dystonia/surfaceFreesurfer/
recon-all -subject s02 -i s02_anatomical.mgh -all -cw256