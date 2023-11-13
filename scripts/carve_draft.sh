#!/bin/bash

carve --refseq GCF_000022065.1 -o ../models/RcH10_draft.xml -u grampos --fbc2 --gapfill MM_xylose,MM_arabinose,MM_glucose,MM_galactose,MM_mannose,MM_cellobiose,DM_cellobiose,DM_cellulose,DM_xyloglucan,DM_cellulose_2,BM_xylan  -i MM_glucose --mediadb ../input/media_db.tsv -v -d --soft ../input/soft_constraints.tsv --hard ../input/hard_constraints.tsv 