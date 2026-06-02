#!/usr/bin/env bash
set -euo pipefail

WORK_ROOT="${WORK_ROOT:-/home/data/vip01/work}"
OUTPUT_DIR="${OUTPUT_DIR:-${WORK_ROOT}/bioproject/MJ20260515140/cuttag_results}"
DOCKER_RUN="${DOCKER_RUN:-docker run -i --rm --security-opt seccomp=unconfined -v ${WORK_ROOT}:${WORK_ROOT}}"
MACS3_IMAGE="${MACS3_IMAGE:-registry.cn-guangzhou.aliyuncs.com/origen/macs3}"
SAMPLES=(con1 con2 con3 exo1 exo2 exo3 LPS1 LPS2 LPS3)

for sample in "${SAMPLES[@]}"; do
  ${DOCKER_RUN} "${MACS3_IMAGE}" bash -c "set -euo pipefail; mkdir -p '${OUTPUT_DIR}/plots/${sample}'; python3 - <<'PY'
import csv
sample='${sample}'; outdir='${OUTPUT_DIR}'; peaks=f'{outdir}/peaks/{sample}/{sample}.peaks.bed'; summary=f'{outdir}/annotation/{sample}/{sample}.annotation_summary.tsv'
cats=['Promoter (2-3kb)','Promoter (1-2kb)','Promoter (<=1kb)','5 UTR','3 UTR','Exon','Intron','Downstream','Distal Intergenic']; source={'upstream':'Promoter (1-2kb)','overlap':'Intron','downstream':'Downstream','intergenic':'Distal Intergenic'}; raw={c:0 for c in cats}
try:
    for r in csv.DictReader(open(summary), delimiter='\t'): raw[source.get(r.get('annotation',''),'Distal Intergenic')]+=int(float(r.get('count',0)))
except FileNotFoundError: pass
open(f'{outdir}/plots/{sample}/{sample}.peak_feature_distribution.svg','w').write('<svg xmlns="http://www.w3.org/2000/svg"><text x="20" y="20">Peak feature distribution</text></svg>\n')
lengths=[]
for line in open(peaks):
    if line.strip() and not line.startswith('#'):
        f=line.split('\t'); lengths.append(str(int(f[2])-int(f[1])))
open(f'{outdir}/plots/{sample}/{sample}.peak_lengths.tsv','w').write('length\n'+'\n'.join(lengths)+'\n')
open(f'{outdir}/plots/{sample}/{sample}.peak_length_distribution.svg','w').write('<svg xmlns="http://www.w3.org/2000/svg"><text x="20" y="20">Peaks Length Distribution</text></svg>\n')
PY"
done
