# spm_toolkit
Python-based toolkit for SPM data manipulation, postprocessing, registrations and visualisation

List of files:
* **spm_data_base.py** - auxiliary constants, variables, routines
* **spm_data_layer.py** - SPMdataLayer class code
* **spm_data_struct.py** - SPMdata class code
* **spm_data_aux.py** - routines working externally with SPMdata class
* **spm_data_inspect.py** - inspect_channels routine and its auxiliary subroutines
* **spm_data_align.py** - routines for alignment and their auxiliary subroutines
* **spm_data_load.py** - load_data routine and its auxiliary subroutines
* **spm_data_cut2d.py** - Cut2D class code
* **spm_data_cut1d.py** - Cut1D class code
* **spm_data_1d.py** - LineData class code
* **registration_base.py** - auxiliary routines for registration
* **registration_class.py** - RegBlock class code
* **registration_compare.py** - CompareRegBlocks class code

**NOTE:** All files assume they are in the same directory. Moreover, this directory contains empty **\_\_init\_\_.py** file so that the directory is a package.
