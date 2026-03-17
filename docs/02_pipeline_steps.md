# Pipeline Workflow and Architecture

The workflow consists of modular scripts located in the `scripts/` directory.

## Pipeline Map (What Each Script Does)

| Step | Script | Function |
| :--- | :--- | :--- |
| **1.** | `01_samtools_bam2fastq.sh` | Converts basecalled `.bam` to `.fastq.gz` |
| **2.** | `02_porechop.sh` | Adapter trimming |
| **3.** | `03_nanofilt.sh` | Length filtering (`<100 bp` removed) |
| **4.** | `03b_read_lengths.sh` | Extrapolates read-length distributions |
| **5.** | `04_nanostat.sh` | Generates QC metrics per sample |
| **6.** | `05_kraken2.sh` | Taxonomic classification on all reads |
| **7.** | `06_listeria_extract.sh` | Extracts target reads (e.g., *Listeria*) and compiles per-sample summaries |
| **8.** | `07_compile_stats.sh` | Compiles read and Target summary tables |
| **9.** | `08_metamdbg.sh`, `08b_myloasm.sh`, `09_metaflye.sh` | Assembly workflows using 3 distinct long-read assemblers |
| **9b.**| `09b_dorado_polish.sh` | Aligns reads back to draft assemblies with `dorado aligner`, then polishes with `dorado polish --bacteria` |
| **10.**| `10_seqkit_fq2fa.sh` | FASTQ to FASTA conversion for reads |
| **11.**| `11_amrfinderplus.sh` | AMR/virulence calls on reads and contigs |
| **12.**| `13_kraken2_contigs.sh` | Taxonomy on assembled contigs |
| **13.**| `14_listeria_contigs.sh` | Extract Target contigs and calculate contig stats |
| **14.**| `15_compile_listeria_overview.sh` | Integrated Target overview compilation |
| **15.**| `16_compile_amr_overview.sh` | AMR overview tables |
| **16.**| `18_assembly_stats.sh` | Assembly summary statistics |
| **17.**| `17_generate_report.sh` | Builds the comprehensive HTML pipeline report |
| **18.**| `20_comparison_report.sh` | Dedicated AS vs N comparative report (optional) |

### Optional Downstream Visualizations/Exports
- `19_reads_report.sh`: Read-focused quick report
- `21_statistical_analysis.py`: statistical testing
- `22_local_plots.py`: Local publication-style figure generation (run this on your laptop!)
- `20_export_tables_to_xlsx.py`, `21_kraken2_to_spreadsheets.py`: Spreadsheet exports

---

## Flag Reference (Core CLI Options Explained)

This section explains the command-line flags uniquely used in our core workflow.

### `samtools` (steps 1 and 9)
- `fastq -@ 4`: `-@` sets worker threads for BAM-to-FASTQ conversion.
- `view -b`: `-b` outputs BAM instead of SAM.

### `porechop` (step 2)
- `-i <file>`: input FASTQ.
- `-o <file>`: output trimmed FASTQ.

### `NanoFilt` (step 3)
- `-l 100`: keep reads with length >= 100 bp.

### `kraken2` (steps 5 and 13)
- `--db <dir>`: Kraken2 database directory.
- `--use-names`: print scientific names in classification output.
- `--report <file>`: write sample-level clade summary report.
- `--output <file>`: write per-read or per-contig classification output (very large!).

### `seqtk` (step 6)
- `subseq <fastq> <id_file>`: extract reads whose IDs are listed in `id_file`.

### `metaMDBG` (step 8)
- `asm`: run assembly mode.
- `--in-ont <fastq>`: ONT reads input.

### `flye` (step 9)
- `--meta`: metagenome assembly mode.
- `--nano-hq <fastq>`: ONT high-quality reads input.

### `minimap2` (step 9)
- `-ax map-ont`: preset for ONT read-to-reference alignment with SAM output.

### `amrfinder` / AMRFinderPlus (step 11)
- `--plus`: include the extended AMRFinderPlus database (AMR + additional marker classes like virulence factors).
- `-n <fasta>`: nucleotide FASTA input.

### `dorado aligner` (step 9b)
- `dorado aligner <draft.fasta> <reads>`: aligns reads (BAM or FASTQ) to a draft assembly reference.
- Output is piped to `samtools sort` and then indexed with `samtools index`.

### `dorado polish` (step 9b)
- `dorado polish <aligned.bam> <draft.fasta>`: runs the polishing model on aligned reads to produce a corrected consensus.
- `--bacteria`: resolves a bacterial-specific polishing model automatically based on the input data type. Use this flag for any bacterial genome assembly.
