# IDMer User Manual
## Usage
iDMer has two modes, `exp` and `denovo` mode.

In `exp` mode, user need supply virus VTPs reliance and restriction gene, EHFs reliance and restriction gene.
In `denovo` mode, iDMer predicted the VTPs based on virus genome information. User need supplies virus genome information and candidate VTPs in fasta format, candidate EHFs reliance and restriction gene.
You can use these two mode by:

        python main.py exp --help

or

        python main.py denovo --help

## Input Files

### Input Files (exp mode)

  VTPs reliance gene  
  VTPs restriction gene  
  EHFs reliance gene  
  EHFs restriction gene  
  for example:  
  python main.py exp -VDN test/VTPs_DN.tsv -VUP teVTPs_UP.tsv -EDN EHFs_DN.tsv  -EUP   EHFs_UP.tsv -o exp
