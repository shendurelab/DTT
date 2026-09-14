
#########################################
### Step-1: creating singlecell.csv sheet
### Python code

data_path = "/net/shendure/vol10/projects/akeith/nobackup/analyses/SR_revisions/ES_multiome/count"
work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape/revision/data_10X_multiome"

import os, sys
import pysam

sample_id = int(sys.argv[1])
print(sample_id)

cell_list = set()
file = open(f"{work_path}/amulet_analysis/ATAC_cell_list_raw.txt")
for line in file:
    l = line.rstrip().split('\t')
    if l[1] == "lane_" + str(sample_id):
        cell_list.add(l[0])
file.close()

barcode_list = set()
bamfile = pysam.AlignmentFile(f"{data_path}/arc_L{sample_id}/outs/atac_possorted_bam.bam", "rb")
for read in bamfile.fetch():
    if read.has_tag("CB"):
        barcode_list.add(read.get_tag("CB"))
bamfile.close()

out = open(f"{work_path}/amulet_analysis/arc_L{sample_id}/singlecell.csv", 'w')
out.write(','.join(['barcode', 'cell_id', 'is__cell_barcode']) + '\n')
for barcode in barcode_list:
    if barcode in cell_list:
        out.write(','.join([barcode, barcode, '1']) + '\n')
    else:
        out.write(','.join([barcode, 'None', '0']) + '\n')
out.close()


#############################
### Step-2: performing AMULET
### Python code

#!/bin/bash
INPUTFILES=(arc_L1 arc_L2 arc_L3 arc_L4)
INPUTFILENAME="${INPUTFILES[$SGE_TASK_ID - 1]}"
echo $INPUTFILENAME

data_path=/net/shendure/vol10/projects/akeith/nobackup/analyses/SR_revisions/ES_multiome/count
work_path=/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape/revision/data_10X_multiome
script_path=/net/shendure/vol10/projects/cxqiu/nobackup/install/AMULET

$script_path/AMULET.sh \
--bambc CB \
--bcidx 0 \
--cellidx 1 \
--iscellidx 2 \
--mapqthresh 30 \
--forcesorted \
$data_path/"$INPUTFILENAME"/outs/atac_possorted_bam.bam \
$work_path/amulet_analysis/"$INPUTFILENAME"/singlecell.csv \
$script_path/mouse_autosomes.txt \
/net/shendure/vol10/projects/cxqiu/nobackup/install/run_AMULET/blacklist_repeats_segdups_rmsk_mm10.bed \
$work_path/amulet_analysis/"$INPUTFILENAME"/ \
$script_path

