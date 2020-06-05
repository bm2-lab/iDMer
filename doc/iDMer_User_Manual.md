# IDMer User Manual
## Usage
iDMer has two modes, `exp` and `denovo` mode.

In `exp` mode, user need supply virus VTPs reliance and restriction genes, EHFs reliance and restriction genes. For example:

        python main.py exp --help
        python main.py exp -VDN  test/VTPs_DN.tsv  -VUP  test/VTPs_UP.tsv  -EDN  test/EHFs_DN.tsv   -EUP  test/EHFs_UP.tsv -o exp
        python main.py exp -VDN  test/VTPs_DN.tsv  -VUP  test/VTPs_UP.tsv  -EDN  test/EHFs_DN.tsv   -EUP  test/EHFs_UP.tsv -o expGAT  --GAT
        
In `denovo` mode, iDMer predicted the VTPs based on virus genome information. User need supply virus genome information and candidate VTPs in fasta format, a config file indicate the VTPs is reliance or restriction gene, candidate EHFs reliance and restriction genes.

        python main.py denovo --help
        python main.py  denovo  -EDN  test/EHFs_DN.tsv  -EUP  test/EHFs_UP.tsv  -virus  test/virus.fa  -host  test/host.fa  
        -config  config.tsv  -o  denovo
        python main.py  denovo  -EDN  test/EHFs_DN.tsv  -EUP  test/EHFs_UP.tsv  -virus  test/virus.fa  -host  test/host.fa   
        -config  config.tsv   -o  denovoGAT  --GAT    

