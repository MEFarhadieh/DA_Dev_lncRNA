### This file contains codes to generate comprehensive lncRNA annotation for alignment, and modified from Singletrome.

# Define input/output directories and parameters
genomesBase=‍‍‍‍~/references
source=$genomesBase/reference_sources
mkdir $source
cd $source

# Download homo sapience reference fasta and GTF files
fasta_url="http://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
fasta_in="${source}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
curl -sS "$fasta_url" | zcat > "$fasta_in"
gtf_url="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.primary_assembly.annotation.gtf.gz"
gtf_in="${source}/gencode.v43.primary_assembly.annotation.gtf"
curl -sS "$gtf_url" | zcat > "$gtf_in"
LncExpDB_gtf_url="https://download.cncb.ac.cn/lncexpdb/0-ReferenceGeneModel/1-GTFFiles/lncRNA_LncBookv1.9_GRCh38.gtf.gz"
LncExpDB_gtf_in="${source}/LncBook_Version2.0_OnlyLnc_hg38.gtf"

cd $source
wget -c $LncExpDB_gtf_url -O - | tar -xz && mv $source/LncBook_Version2.0_onlylnc.gtf $LncExpDB_gtf_in

# Remove genes with invalid exons order in trancript
correctedLncGtf="${source}/$(basename "$LncExpDB_gtf_in")_corrected.gtf"
cat $LncExpDB_gtf_in | grep -v 'HSALNG0056858\|HSALNG0059740\|HSALNG0078365\|HSALNG0092690\|HSALNG0093062\|HSALNG0089130\|HSALNG0089954\|HSALNG0095105'> $correctedLncGtf

# Modify sequence headers in the Ensembl FASTA
fasta_modified="$source/$(basename "$fasta_in").modified"
cat "$fasta_in" \
    | sed -E 's/^>(\S+).*/>\1 \1/' \
    | sed -E 's/^>([0-9]+|[XY]) />chr\1 /' \
    | sed -E 's/^>MT />chrM /' \
    > "$fasta_modified"
gtf_modified="$source/$(basename "$gtf_in").modified"
ID="(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"
cat "$gtf_in" \
    | sed -E 's/gene_id "'"$ID"'";/gene_id "\1"; gene_version "\3";/' \
    | sed -E 's/transcript_id "'"$ID"'";/transcript_id "\1"; transcript_version "\3";/' \
    | sed -E 's/exon_id "'"$ID"'";/exon_id "\1"; exon_version "\3";/' \
    > "$gtf_modified"

# Define string patterns for GTF tags
BIOTYPE_PATTERN="(protein_coding)"
GENE_PATTERN="gene_type \"${BIOTYPE_PATTERN}\""
TX_PATTERN="transcript_type \"${BIOTYPE_PATTERN}\""
READTHROUGH_PATTERN="tag \"readthrough_transcript\""
PAR_PATTERN="tag \"PAR\""
cat "$gtf_modified" \
    | awk '$3 == "transcript"' \
    | grep -E "$GENE_PATTERN" \
    | grep -E "$TX_PATTERN" \
    | grep -Ev "$READTHROUGH_PATTERN" \
    | grep -Ev "$PAR_PATTERN" \
    | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
    | sort \
    | uniq \
    > "${source}/gene_allowlist_protein_coding"
gtf_protein_coding="${source}/$(basename "$gtf_in").protein_coding.gtf"
grep -E "^#" "$gtf_modified" > "$gtf_protein_coding"
grep -Ff "${source}/gene_allowlist_protein_coding" "$gtf_modified" >> "$gtf_protein_coding"

# Genrate a GTF for Gencode lncRNAs only for overlapping with Gencode protein coding
BIOTYPE_PATTERNLncRNA="(lncRNA)"
GENE_PATTERNLncRNA="gene_type \"${BIOTYPE_PATTERNLncRNA}\""
cat "$gtf_modified" | awk '$3 == "transcript"' | grep -E "$GENE_PATTERNLncRNA" | sed -E 's/.*(gene_id "[^"]+").*/\1/' | sort | uniq > "${source}/gene_allowlist_gencode_lncRNA"
gtf_geneCodeLncRNAPath="${source}/$(basename "$gtf_in").gencode_lncRNA.gtf"
grep -E "^#" "$gtf_modified" > "$gtf_geneCodeLncRNAPath"
grep -Ff "${source}/gene_allowlist_gencode_lncRNA" "$gtf_modified" >> "$gtf_geneCodeLncRNAPath"
