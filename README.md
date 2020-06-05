# iDMer: an integrative Data and Mechanism-driven epidemic response system to identify compound interventions for sudden virus outbreak

iDMer is presented to identify compound interventions for any virus as long as the viral genome is sequenced. It consists of three mechanism-oriented compound identification modules, i.e. the virus-host interaction-oriented module, the autophagy-oriented module, and the cytokine release syndrome-oriented module. In addition, the evaluation of the predicted compound toxicities and the identification of compound combinations for virus treatment with clear mechanisms are also incorporated into iDMer as a one-stop integrative platform.

#### Authors:
Zhiting Wei, Yuli Gao and Fangliangzi Meng

#### Citation:
iDMer: an integrative Data and Mechanism-driven epidemic response system to identify compound interventions for sudden virus outbreak

## Dependencies

#### Required Software:
* [dgl](https://www.dgl.ai/)
* [pytorch](https://pytorch.org/)
* [scikit-learn](https://scikit-learn.org/stable/index.html)
* [miniconda2](https://docs.conda.io/en/latest/miniconda.html)
* [miniconda3](https://docs.conda.io/en/latest/miniconda.html)
* [cmapPy](https://clue.io/cmapPy/index.html)
* [deepDDI](https://bitbucket.org/kaistsystemsbiology/deepddi/src/master/)
* [ProTox-II](http://tox.charite.de/protox_II)
* [HVPPI](http://zzdlab.com/hvppi/)
* [biopython](https://biopython.org/)

## Installation
#### Install via docker, highly recommended
Docker image of iDMer is available at https://hub.docker.com/r/bm2lab/idmer/.
if you have docker installed, you can pull the image like so:

        docker pull bm2lab/idmer

#### Install from source, not recommended

        git clone https://github.com/bm2-lab/iDMer.git

## Usage
iDMer has two modes, `exp` and `denovo` mode.

In `exp` mode, user need supply virus VTPs reliance and restriction gene, EHFs reliance and restriction gene.
In `denovo` mode, iDMer predicted the VTPs based on virus genome information. User need supplies virus genome information and candidate VTPs in fasta format, candidate EHFs reliance and restriction gene.
You can use these two mode by:

        python main.py exp --help

or

        python main.py denovo --help

## User Manual
For detailed information about usage, input and output files, test examples and data preparation please refer to the [iDMer User Manual](/doc/iDMer_User_Manual.md)

## Contact
Zhiting Wei 1632738@tongji.edu.cn  
Qi Liu qiliu@tongji.edu.cn  
Tongji University, Shanghai, China
