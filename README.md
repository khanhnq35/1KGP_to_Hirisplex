# 🧬 HiRisPlex-S Genotype Extraction Pipeline (1KGP Phase 3)

## 📘 Overview
This project automates the extraction of **41 HiRisPlex-S SNPs** from the **1000 Genomes Project Phase 3** VCFs  
and converts them into a **per-sample genotype matrix (REF/ALT-aware)** suitable for downstream eye-/hair-/skin-color prediction.

```
├── extract_1kgp_make_matrix.sh
├── make_hirisplex_by_sample_refaware.py
├── hirisplexs.csv                ← official HiRisPlex-S SNP panel
├── 1kgp_out/
│   ├── rs41_ids.txt
│   ├── rs41_GT_matrix_with_refalt.tsv
│   ├── samples_header.tsv
│   └── hirisplex_by_sample_refaware.csv
└── Result_Final.csv              ← your later analysis output
```

---

## ⚙️ Environment
- **OS** : macOS (bash 3.2 +)
- **Tools required**
  - [`bcftools`](http://samtools.github.io/bcftools/)
  - `curl`, `awk`, `sort`, `head`, `paste`
  - Python ≥ 3.8 (with standard library only)

> 💡 Install bcftools via Homebrew  
> ```bash
> brew install bcftools
> ```

---

## 🧩 Step 1 – Extract 41 SNPs from 1KGP

Run once to download and merge all 41 SNPs from chromosomes 1–22 + X.

```bash
chmod +x extract_1kgp_make_matrix.sh
./extract_1kgp_make_matrix.sh
```

### 📤 Output
| File | Description |
|------|--------------|
| `1kgp_out/rs41_ids.txt` | List of 41 rsIDs (HiRisPlex-S panel) |
| `1kgp_out/samples_header.tsv` | Sample IDs from 1000 Genomes |
| `1kgp_out/rs41_GT_matrix_with_refalt.tsv` | Combined genotype matrix (`CHROM POS ID REF ALT GTs…`) |

> The script safely resumes downloads, uses `bcftools norm -m -both` to split multi-allelic records,  
> and logs progress for each chromosome.

---

## 🧠 Step 2 – Convert Matrix → HiRisPlex-S CSV (REF-Aware)

```bash
python3 make_hirisplex_by_sample_refaware.py
```

### 📥 Input
- `1kgp_out/rs41_GT_matrix_with_refalt.tsv`  
- `hirisplexs.csv` (the official panel listing SNP ↔ Allele)

### 📤 Output
| File | Description |
|------|--------------|
| `1kgp_out/hirisplex_by_sample_refaware.csv` | Final CSV matrix: `sampleid, rsID_Allele, …` |
| Each value = `0 / 1 / 2 / NA` (count of target allele per sample) |

The Python script:
- Reads the panel → maps rsID ↔ target allele  
- Handles REF/ALT/complement and simple INDEL logic  
- Outputs per-sample rows in UTF-8 CSV  
- Logs missing or non-matching alleles (`⚠️ rs1805008_T toàn NA…`)  

---

## 🧪 Example Usage

Preview 3 lines of the merged matrix:
```bash
head -n 3 1kgp_out/rs41_GT_matrix_with_refalt.tsv
```

Preview 3 lines of the final CSV:
```bash
head -n 3 1kgp_out/hirisplex_by_sample_refaware.csv
```

---

## 🧼 Notes
- `bcftools norm -m -both` ensures multi-allelic records (e.g., ALT =A,C,T) are split and GT re-coded.  
- Missing alleles not present in 1KGP are counted as `0` or `NA`.  
- You can rerun Python script without re-extracting VCFs.

---

## 🧾 Citation
HiRisPlex-S System — Walsh S et al., *Forensic Science International: Genetics* (2018).  
1000 Genomes Project — *Nature* 526, 68–74 (2015).

---

## 👨‍💻 Author
Nguyễn Quốc Khánh (20235118) — SOICT @ HUST  
For questions: <your email or GitHub link here>
