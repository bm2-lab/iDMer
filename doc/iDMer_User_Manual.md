# iDMer User Manual
## Usage
iDMer has two modes, `exp` and `denovo` mode.
## exp mode
In `exp` mode, user needs to provide virus VTPs reliance and restriction genes, EHFs reliance and restriction genes. As described in our manuscript, if the test virus has been extensively investigated, VTPs and EHFs can be collected from the Viruses-STRING database and literatures. Subsequently, please classify the curated VTPs and EHFs of the test virus as reliance or restriction genes accroding to their functions described in KEGG and Reactome Pathway Database. Finally, run iDMer exp mode as follows:

        cd /home/test
        python /home/main.py exp --help
        python /home/main.py exp -VDN  VTPs_DN.tsv  -VUP  VTPs_UP.tsv  -EDN  EHFs_DN.tsv   -EUP  EHFs_UP.tsv -output exp
        or
        python /home/main.py exp -VDN  VTPs_DN.tsv  -VUP  VTPs_UP.tsv  -EDN  EHFs_DN.tsv   -EUP  EHFs_UP.tsv -output expGAT  --GAT  

### Output Files
down_proteins: valid and deduplicated virus reliance genes  
up_proteins: valid and deduplicated virus restriction genes  
up_down_protein_GAT.csv: ranked host factors output by the GAT model   
CMap_tox.tsv: ranked compound interventions with toxicity annotation   


## denovo mode
In `denovo` mode, (1) As for the test virus' VTPs, iDMer applies HVPPI to predict them based on the test virus genome information and the candidate VTPs that may interact with the test virus. The set of union VTPs from viruses that belongs to the same family of the test virus are taken as the candidate VTPs. Subsequently, please classify the predicted VTPs as reliance or restriction genes accroding to their functions described in KEGG and Reactome Pathway Database and make a config file to indicate whether the predicted VTPs are reliance or restriction genes; (2) As for the test virus' EHFs, the set of overlapping EHFs from viruses that belongs to the same family of the test virus are taken as the test virus' EHFs. Classify the test virus' EHFs as reliance or restriction genes accroding to their functions described in KEGG and Reactome Pathway Database. Finally, run iDMer denovo mode as follows:

        cd /home/test
        python /home/main.py  denovo --help
        python /home/main.py  denovo  -EDN  EHFs_DN.tsv  -EUP  EHFs_UP.tsv  -virus  virus.fa  -host  host.fa  -config  config.tsv  -output  denovo
        or
        python /home/main.py  denovo  -EDN  EHFs_DN.tsv  -EUP  EHFs_UP.tsv  -virus  virus.fa  -host  host.fa  -config  config.tsv  -output denovoGAT --GAT    

### Output Files
VTPs_DN.tsv: predicted VTPs reliance genes  
VTPs_UP.tsv: predicted VTPs restriction genes   
down_proteins: valid and deduplicated virus reliance genes     
up_proteins: valid and deduplicated virus restriction genes     
up_down_protein_GAT.csv: ranked host factors output by the GAT model  
CMap_tox.tsv: ranked compound interventions with toxicity annotation  

### Column explanation
The output CMap_tox.tsv consists of the following columns:

| Column Name           | Description |
| -----------           | ----------- |
| pert_id               | id in the CMap Touchstone |
| TAG                   | Compound CMap score |
| pert_iname            | Compound name |
| Toxicity              | Whether the identified compound is toxic or not|

## Compound combination identification
Compound identified in the above step can be combined with compound identified in the CRS-oriented module. We need to activate the deepDDI conda environment firstly. For example,

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
| pert_iname            | Compound name |
| ddi_type              | Interaction type between two compounds |
| sentence              | Interaction type description
| antagonism            | Whether a compound combination is antagonistic or not |
| combine score         | mean CMap score for two compounds









