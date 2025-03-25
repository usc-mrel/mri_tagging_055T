# MRI Tagging at 0.55T

This is a collection of pulse sequences developed in PyPulseq, and reconstruction algorithms written in MATLAB.

## Dataset
Dataset is located [here](https://zenodo.org/records/15079693).

## Pulse Sequence Design
Pulse Sequences are developed using Pypulseq and are located as a submodule to another repository in sequences/rtspiral_pypulseq. To generate pulse sequences, please install all dependencies using the README located in [rtspiral_pypulseq](https://github.com/usc-mrel/rtspiral_pypulseq/tree/1721c6efdcb8dc940a0cfe7b1fd642068fe994b8).

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
