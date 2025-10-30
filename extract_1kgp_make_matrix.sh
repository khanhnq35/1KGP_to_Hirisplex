#!/usr/bin/env bash
# ============================================================
#  Script: extract_1kgp_make_matrix.sh
#  Purpose: Extract 41 HiRisPlex-S SNPs from 1KGP VCFs
#  Output : ~/Documents/2025.1/Bio_Labs/1kgp_out/rs41_GT_matrix_with_refalt.tsv
#  Compatible: macOS bash 3.2
#  ------------------------------------------------------------
#  Usage:
#      chmod +x extract_1kgp_make_matrix.sh
#      ./extract_1kgp_make_matrix.sh
# ============================================================

set -euo pipefail
trap 'echo "❌ Error at line $LINENO. Exit code $?."' ERR

# ==== CONFIG ====
BASE="${HOME}/Documents/2025.1/Bio_Labs"
OUTDIR="${BASE}/1kgp_out"
TMPDIR="${BASE}/1kgp_tmp"
mkdir -p "${OUTDIR}" "${TMPDIR}"

MIRROR1="https://ftp-trace.ncbi.nlm.nih.gov/1000genomes/ftp/release/20130502"
MIRROR2="https://ftp.ebi.ac.uk/pub/databases/1000genomes/ftp/release/20130502"

OUTFILE="${OUTDIR}/rs41_GT_matrix_with_refalt.tsv"
RSFILE="${OUTDIR}/rs41_ids.txt"
HDRFILE="${OUTDIR}/samples_header.tsv"

# ==== STEP 1: Create rs41_ids.txt if missing ====
if [ ! -s "${RSFILE}" ]; then
  echo "[+] Creating ${RSFILE}"
  cat > "${RSFILE}" <<'RS41'
rs312262906
rs11547464
rs885479
rs1805008
rs1805005
rs1805006
rs1805007
rs1805009
rs201326893
rs2228479
rs1110400
rs28777
rs16891982
rs12821256
rs4959270
rs12203592
rs1042602
rs1800407
rs2402130
rs12913832
rs2378249
rs12896399
rs1393350
rs683
rs3114908
rs1800414
rs10756819
rs2238289
rs17128291
rs6497292
rs1129038
rs1667394
rs1126809
rs1470608
rs1426654
rs6119471
rs1545397
rs6059655
rs12441727
rs3212355
rs8051733
RS41
else
  echo "[=] Using existing ${RSFILE}"
fi

# ==== STEP 2: Determine sample header ====
echo "[+] Fetching sample list for header..."
VCF1="${TMPDIR}/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
VCFURL1="${MIRROR1}/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

if [ ! -s "${VCF1}" ]; then
  echo "  ↪ Downloading small header from NCBI..."
  if ! curl -fsI "${VCFURL1}" >/dev/null 2>&1; then
    VCFURL1="${MIRROR2}/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
  fi
  bcftools query -l "${VCFURL1}" > "${HDRFILE}" 2>/dev/null || {
    echo "❌ Failed to get sample list from mirrors."
    exit 1
  }
else
  bcftools query -l "${VCF1}" > "${HDRFILE}"
fi

if [ ! -s "${HDRFILE}" ]; then
  echo "❌ Cannot obtain sample list. Abort."
  exit 1
fi

# Build header line
header="CHROM\tPOS\tID\tREF\tALT"
while read s; do header="${header}\t${s}"; done < "${HDRFILE}"
echo -e "${header}" > "${OUTFILE}"
echo "[+] Header -> ${HDRFILE}"
echo "[+] Output initialized: ${OUTFILE}"

# ==== STEP 3: Define function for suffix ====
get_suffix() {
  local chr="$1"
  if [ "${chr}" = "X" ]; then
    echo "v1b"
  else
    echo "v5a"
  fi
}

# ==== STEP 4: Main extraction loop ====
for CHR in $(eval echo {1..22}) X; do
  echo "==> Processing chr${CHR}"
  SUFFIX=$(get_suffix "${CHR}")
  FN="ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_${SUFFIX}.20130502.genotypes.vcf.gz"
  LOCAL="${TMPDIR}/${FN}"
  URL1="${MIRROR1}/${FN}"
  URL2="${MIRROR2}/${FN}"
  TMP="${OUTDIR}/_chr${CHR}_rs41.tmp.tsv"

  SRC=""
  if [ -s "${LOCAL}" ]; then
    SRC="${LOCAL}"
    echo "  ✅ Using local file: ${SRC}"
  else
    echo "  ↪ Attempting download..."
    if curl -fsI "${URL1}" >/dev/null 2>&1; then
      SRC="${URL1}"
    elif curl -fsI "${URL2}" >/dev/null 2>&1; then
      SRC="${URL2}"
    else
      echo "  ⚠️  Both mirrors failed for chr${CHR}, skipping."
      continue
    fi
  fi

  echo "  ▶ Extracting SNPs..."
  if bcftools view -i "ID=@${RSFILE}" "${SRC}" -Ou 2>/dev/null | \
     bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' > "${TMP}" 2>/dev/null; then
     LINES=$(wc -l < "${TMP}" | awk '{print $1}')
     echo "  ✅ Extracted ${LINES} lines."
     cat "${TMP}" >> "${OUTFILE}"
  else
     echo "  ❌ Extraction failed for chr${CHR}"
     continue
  fi
  rm -f "${TMP}"
  echo "[DONE] chr${CHR}"
done

# ==== STEP 5: Show preview ====
echo
echo "✅ All done!"
echo "[+] Output file: ${OUTFILE}"
echo "[+] Preview (first 3 lines):"
head -n 3 "${OUTFILE}"