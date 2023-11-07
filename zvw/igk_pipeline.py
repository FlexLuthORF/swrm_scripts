import os
import subprocess
import sys

# Define functions

def align_with_minimap2(fasta, prefix, reffn, threads):
    # Align sequences using minimap2 and perform sorting and indexing
    subprocess.run(f'minimap2 -t {threads} -L -a {reffn} {fasta} > {prefix}.sam', shell=True)
    subprocess.run(f'samtools view -Sbh {prefix}.sam > {prefix}.bam', shell=True)
    subprocess.run(f'samtools sort -@ {threads} {prefix}.bam -o {prefix}.sorted.bam', shell=True)
    subprocess.run(f'samtools index {prefix}.sorted.bam', shell=True)
    subprocess.run(f'rm -f {prefix}.sam', shell=True)
    subprocess.run(f'rm -f {prefix}.bam', shell=True)

def align_with_minimap2_asm20(fasta, prefix, reffn, threads):
    # Align sequences using minimap2 with asm20 preset and perform sorting and indexing
    subprocess.run(f'minimap2 -x asm20 -t {threads} -L -a {reffn} {fasta} > {prefix}.sam', shell=True)
    subprocess.run(f'samtools view -Sbh {prefix}.sam > {prefix}.bam', shell=True)
    subprocess.run(f'samtools sort -@ {threads} {prefix}.bam -o {prefix}.sorted.bam', shell=True)
    subprocess.run(f'samtools index {prefix}.sorted.bam', shell=True)
    subprocess.run(f'rm -f {prefix}.sam', shell=True)
    subprocess.run(f'rm -f {prefix}.bam', shell=True)

def remove_softclips():
    # Remove soft clips from alignments
    for iter in ['1', '2']:
        os.makedirs(f'{outdir}/break_at_soft_clip/{iter}', exist_ok=True)
        for i in ['1', '2']:
            bam = f'{outdir}/hifiasm/asm.bp.hap{i}.p_ctg_to_ref.sorted.bam' if iter == '1' else f'{outdir}/break_at_soft_clip/1/{i}_hifi_asm_to_ref.sorted.bam'
            if not os.path.exists(f'{outdir}/break_at_soft_clip/{iter}/{i}_hifi_asm_to_ref.sorted.bam'):
                subprocess.run(f'python /home/egenge01/bioinformatics/bioinformatics-common/python/extract_soft_clip_seq.py {bam} > {outdir}/break_at_soft_clip/{iter}/{i}_hifi_asm.fasta', shell=True)
                subprocess.run(f'samtools faidx {outdir}/break_at_soft_clip/{iter}/{i}_hifi_asm.fasta', shell=True)
                align_with_minimap2(
                    f'{outdir}/break_at_soft_clip/{iter}/{i}_hifi_asm.fasta',
                    f'{outdir}/break_at_soft_clip/{iter}/{i}_hifi_asm_to_ref',
                    reffn,
                    threads
                )
            if not os.path.exists(f'{outdir}/break_at_soft_clip/{iter}/{i}_asm20_hifi_asm_to_ref.sorted.bam'):
                align_with_minimap2_asm20(
                    f'{outdir}/break_at_soft_clip/{iter}/{i}_hifi_asm.fasta',
                    f'{outdir}/break_at_soft_clip/{iter}/{i}_asm20_hifi_asm_to_ref',
                    reffn,
                    threads
                )

# Main part of the script

outdir = sys.argv[1]
ccs = sys.argv[2]
threads = 10

reffn = "/home/egenge01/projects/IGK/data/make_franken/franken_ref/reference.fasta"

if not os.path.exists(f'{outdir}/reads.fasta.fai'):
    subprocess.run(f'samtools view {ccs} | awk \'{{ print ">"$1"\\n"$10 }}\' > {outdir}/reads.fasta', shell=True)
    subprocess.run(f'samtools faidx {outdir}/reads.fasta', shell=True)

if not os.path.exists(f'{outdir}/hifiasm/asm.bp.hap2.p_ctg.fasta.fai'):
    os.makedirs(f'{outdir}/hifiasm', exist_ok=True)
    subprocess.run(f'hifiasm -o {outdir}/hifiasm/asm -t {threads} {outdir}/reads.fasta', shell=True)
    subprocess.run(f'gfatools gfa2fa {outdir}/hifiasm/asm.bp.p_utg.gfa > {outdir}/hifiasm/asm.bp.p_utg.fasta', shell=True)
    for i in ['1', '2']:
        subprocess.run(f'gfatools gfa2fa {outdir}/hifiasm/asm.bp.hap{i}.p_ctg.gfa > {outdir}/hifiasm/asm.bp.hap{i}.p_ctg.fasta', shell=True)
        subprocess.run(f'samtools faidx {outdir}/hifiasm/asm.bp.hap{i}.p_ctg.fasta', shell=True)

for i in ['1', '2']:
    fn = f'asm.bp.hap{i}.p_ctg'
    if not os.path.exists(f'{outdir}/hifiasm/{fn}_to_ref.sorted.bam.bai'):
        align_with_minimap2_asm20(
            f'{outdir}/hifiasm/asm.bp.hap{i}.p_ctg.fasta',
            f'{outdir}/hifiasm/{fn}_to_ref',
            reffn,
            threads
        )

for i in ['p']:
    fn = f'asm.bp.{i}_utg'
    if not os.path.exists(f'{outdir}/hifiasm/{fn}_to_ref.sorted.bam.bai'):
        align_with_minimap2_asm20(
            f'{outdir}/hifiasm/{fn}.fasta',
            f'{outdir}/hifiasm/{fn}_to_ref',
            reffn,
            threads
        )

# Call the function to remove soft clips
remove_softclips()
