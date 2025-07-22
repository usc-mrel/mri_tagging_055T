# MRI Tagging at 0.55T

This is a collection of pulse sequences developed in PyPulseq, and reconstruction algorithms written in MATLAB.

## Dataset
Dataset is located [here](https://zenodo.org/records/15079693). Please place the data in the `data` directory in the root of this folder.

## Pulse Sequence Design
Pulse Sequences are developed using Pypulseq and are located as a submodule to another repository in sequences/rtspiral_pypulseq. To generate pulse sequences, please install all dependencies using the README located in [rtspiral_pypulseq](https://github.com/usc-mrel/rtspiral_pypulseq/tree/1721c6efdcb8dc940a0cfe7b1fd642068fe994b8). We tested all sequence development code using the python version 3.11.8.

To run pulse sequences, please `cd sequences/rtspiral_pypulseq`

To generate the real-time speech sequence, please run 
```
python write_rtspiral.py -c protocols/speech_bssfp_dynamic_tagging_055T.toml
```

To generate real-time cardiac sequence, please run 
```
python write_rtspiral_cine.py -c protocols/cardiac_bssfp_cine_055T.toml
``` 

Expected output: sequences will be generated in the `out_seq` folder.

## Reconstruction
Reconstruction depends on the `usc_dynamic_reconstruction` toolbox, `ismrmrd`, and `ismrm_sunrise_matlab`, which are all included as submodules in this repository. 
In order to use the `usc_dynamic_reconstruction` toolbox, please follow it's README. To use it, the [Michigan Image Reconstruction Toolbox (MIRT)](https://github.com/JeffFessler/mirt) is requied as a dependency.

### Usage:
By default, the use of the GPU is set to off, but if you have a GPU, turning it on will greatly increase reconstruction speed. Please change the flag:
```
USE_GPU=1 % (line 19 on either recon script)
```
to enable.

To reconstruct speech data, you can run the following in your MATLAB terminal after downloading the data: <br>
```
cd reconstruction/speech
recon_STCR_2D
```

Expected output: reconstruction in `reconstruction/speech/recon/*.mat` inside image_stcr variable.

To reconstruct speech data, you can run the following in your MATLAB terminal: <br>
```
cd reconstruction/cardiac
recon_cine_2d
```
Expected output: reconstruction in `reconstruction/speech/recon/*.mat` inside image_stcr variable.
  
## Directory Hierarchy for Raw Data
The code expects raw data in the following hierarchy:

    DATA_ROOT\
        |- data_folder\
            |- SEQUENCEHASH.mat
            |- raw\
                |- h5\
                    |- raw_file.h5
                |- noise\
                    |- noise_raw_file.h5

`SEQUENCEHASH.mat` is the metadata file generated during sequence design. Refer to [rtspiral_pypulseq](https://github.com/usc-mrel/rtspiral_pypulseq) for the details.


SEQUENCEHASH.mat is the metadata file generated during sequence design, and should be put in `reconstruction/trajectories/`. Refer to rtspiral_pypulseq for the details.


## Attributions

This code includes elements from various sources, and a best effort has been made to attribute all contributions accurately. If you identify any missing or incorrect attributions, please contact me (prakashk AT usc DOT edu) so that we can correct them. We appreciate your understanding and cooperation in ensuring proper attribution for all utilized code.

This project includes code from the following open source libraries: 
- Pypulseq [AGPL](https://github.com/imr-framework/pypulseq/blob/master/LICENSE)
- ISMRM sunrise toolbox [no license](https://github.com/hansenms/ismrm_sunrise_matlab)
- Michigan IRT [MIT](https://github.com/JeffFessler/mirt/blob/main/LICENSE)
- usc_dynamic_reconstruction [MIT](https://github.com/usc-mrel/usc_dynamic_reconstruction/blob/main/LICENSE)

Copyright (c) 2025 Prakash Kumar, Magnetic Resonance Engineering Laboratory. 
[MIT license](https://github.com/usc-mrel/mri_tagging_055T/blob/main/LICENSE).


