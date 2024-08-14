
BEDPE=/home/database/data/InferX/data/frontal_cortex/scATAC/rds/ciceroFrame_conns_uniq.bedpe

COOL=/home/database/data/InferX/data/frontal_cortex/scATAC/scHIC/Astro_all_brain.txt_1kb_contacts.cool
OUT=$BEDPE\.Astro.count
nohup /home/toolkit/local/bin/python3 coolLoopCount.py $COOL $BEDPE $OUT &

COOL=/home/database/data/InferX/data/frontal_cortex/scATAC/scHIC/MG_all_brain.txt_1kb_contacts.cool
OUT=$BEDPE\.MG.count
nohup /home/toolkit/local/bin/python3 coolLoopCount.py $COOL $BEDPE $OUT &

COOL=/home/database/data/InferX/data/frontal_cortex/scATAC/scHIC/Neuron.Ex_all_brain.txt_1kb_contacts.cool
OUT=$BEDPE\.Neuron.Ex.count
nohup /home/toolkit/local/bin/python3 coolLoopCount.py $COOL $BEDPE $OUT &

COOL=/home/database/data/InferX/data/frontal_cortex/scATAC/scHIC/Neuron.In_all_brain.txt_1kb_contacts.cool
OUT=$BEDPE\.Neuron.In.count
nohup /home/toolkit/local/bin/python3 coolLoopCount.py $COOL $BEDPE $OUT &

COOL=/home/database/data/InferX/data/frontal_cortex/scATAC/scHIC/ODC_all_brain.txt_1kb_contacts.cool
OUT=$BEDPE\.ODC.count
nohup /home/toolkit/local/bin/python3 coolLoopCount.py $COOL $BEDPE $OUT &

COOL=/home/database/data/InferX/data/frontal_cortex/scATAC/scHIC/OPC_all_brain.txt_1kb_contacts.cool
OUT=$BEDPE\.OPC.count
nohup /home/toolkit/local/bin/python3 coolLoopCount.py $COOL $BEDPE $OUT &
