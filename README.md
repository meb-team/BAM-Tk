# bamtk

## Description:

BAM-Tk is a software toolkit for dealing with Binary Alignment Map (BAM) files. 
It allows : 
1. abundance matrix construction 
2. *futures functionalities*


## Running BAM-Tk 

### Python libraries
BAM-Tk is designed for Python >=3.6 and requires the following libraries, which will be automatically installed:
* biolib >=0.1.6: common tasks in bioinformatic
* NumPy >=1.9.0: scientific computing with Python.

### Third-party software
BAM-Tk makes use of the following 3rd party dependencies and assumes they are on your system path:
 * [samtools](https://github.com/samtools/samtools) >= 0.1.19: Li H., et al. 2009 The Sequence alignment/map (SAM) format and SAMtools Bioinformatics, 25, 2078-9.

Please cite these tools if you use BAM-Tk in your work.

### Installation

Once the third-party dependencies have been installed, BAM-Tk can be installed using pip:

```
python -m pip install git+https://github.com/meb-team/BAM-Tk.git
```


### Need help ? 

1. Access the help menu (`bamtk -h`)
2. Run an example (`bamtk mm_features tests/data/toy.fa.fai tests/bamlist.txt tests/results`)

## Bugs

* Submit problems or requests here: https://github.com/meb-team/BAM-Tk/issues/

## Citation

Written by Corentin Hochart (corentin.hochart@uca.fr), UMR CNRSS 6023 Laboratoire Genome et Environement (LMGE) as part of the [ANR Eureka](https://anr.fr/Projet-ANR-14-CE02-0004) project. 
Released under the terms of the GNU General Public License v3. bamtk version v0.1.1.

