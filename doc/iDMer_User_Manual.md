# IDMer User Manual
## Usage
iDMer has two modes, `exp` and `denovo` mode.
## exp mode
In `exp` mode, user need supply virus VTPs reliance and restriction genes, EHFs reliance and restriction genes. For example:

        cd /home/test
        python /home/main.py exp --help
        python /home/main.py exp -VDN  VTPs_DN.tsv  -VUP  VTPs_UP.tsv  -EDN  EHFs_DN.tsv   -EUP  EHFs_UP.tsv -output exp
        or
        python /home/main.py exp -VDN  VTPs_DN.tsv  -VUP  VTPs_UP.tsv  -EDN  EHFs_DN.tsv   -EUP  EHFs_UP.tsv -output expGAT  --GAT  

### Output Files
down_proteins: valid and deduplicated virus reliance gene  
up_proteins: valid and deduplicated virus restriction gene  
up_down_protein_GAT.csv: ranked host factors output by the GAT model   
CMap_tox.tsv: ranked compound interventions with toxicity annotation   


## denovo mode
In `denovo` mode, iDMer predicted the VTPs based on virus genome information. User need supply virus genome information and candidate VTPs in fasta format, a config file indicate the VTPs is reliance or restriction gene, candidate EHFs reliance and restriction genes. For examples:

        cd /home/test
        python /home/main.py  denovo --help
        python /home/main.py  denovo  -EDN  EHFs_DN.tsv  -EUP  EHFs_UP.tsv  -virus  virus.fa  -host  host.fa  -config  config.tsv  -output  denovo
        or
        python /home/main.py  denovo  -EDN  EHFs_DN.tsv  -EUP  EHFs_UP.tsv  -virus  virus.fa  -host  host.fa  -config  config.tsv  -output denovoGAT --GAT    

### Output Files
VTPs_DN.tsv: predicted VTPs reliance gene  
VTPs_UP.tsv: predicted VTs restriction gene    
down_proteins: valid and deduplicated virus reliance gene     
up_proteins: valid and deduplicated virus restriction gene     
up_down_protein_GAT.csv: ranked host factors output by the GAT model  
CMap_tox.tsv: ranked compound interventions with toxicity annotation  

### Column explanation
The output CMap_tox.tsv consists of the following columns:

| Column Name           | Description |
| -----------           | ----------- |
| pert_id               | id in the CMap Touchstone |
| TAG                   | Compound CMap score |
| pert_iname            | Compound common name |
| Toxicity              | Whether the identified compound is toxic or not|

## Compound combination identification
compound identified in the above step can be combine with compound identified in the CRS-oriented module. First, we need to activate the deepDDI conda environment. For example,

        conda activate deepDDI
        cd  /home/test
        python  /home/comb.py --help
        python  /home/comb.py  exp/CMap_tox.tsv exp
        conda deactivate
### Output Files
compoundCombination.tsv

### Column explanation
The output compoundCombination.tsv consists of the following columns:

| Column Name           | Description |
| -----------           | ----------- |
| pert_id               | id in the CMap Touchstone |
| pert_iname            | Compound common name |
| ddi_type              | Interaction type between the two compounds |
| sentence              | Interaction type description
| antangonism           | Whether a compound combination is antangonism or not |
| combine score         | Compound combination average CMap score









