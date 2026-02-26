# Installation and Setup

This pipeline is designed for Linux/macOS with conda-compatible package management (Conda or Mamba). For Windows, we recommend using WSL2.

## 1) Install Mamba (Recommended for Speed)
Mamba is a drop-in replacement for Conda that resolves environments significantly faster.

**If you don't have Mamba installed yet:**
```bash
# Example for installing Miniforge (which includes mamba) on Linux:
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```

## 2) Set up Bioconda channels (one time)
```bash
mamba config --add channels bioconda
mamba config --add channels conda-forge
mamba config --set channel_priority strict
```

## 3) Create an environment with all required tools
The pipeline uses one consolidated environment for all steps.

```bash
mamba create -n listeria_as \
  python=3.10 \
  samtools porechop nanofilt nanostat kraken2 seqtk seqkit \
  flye metamdbg myloasm minimap2 racon ncbi-amrfinderplus \
  pandas numpy scipy matplotlib \
  -c conda-forge -c bioconda --strict-channel-priority
```

Then activate it:
```bash
mamba activate listeria_as
```
*(Note: If the scripts on your cluster specifically look for `conda activate tim`, you can either name your environment `tim` or update the `# Activate conda environment` section at the top of the `scripts/*.sh` files.)*

### Quick check that tools are available:
```bash
samtools --version
kraken2 --version
flye --version
amrfinder --version
```

## 4) Databases: Where to get them and how to configure

### Kraken2 database (required for steps 5 and 13)
Current scripts use a `KRAKEN2_DB` variable, which you must set in:
- `scripts/05_kraken2.sh`
- `scripts/13_kraken2_contigs.sh`

**Options to get the database:**
1. **Prebuilt (Easiest):** Download the standard or PlusPF index from the [Ben Langmead AWS index collection](https://benlangmead.github.io/aws-indexes/k2).
2. **Build from NCBI/RefSeq:**
```bash
kraken2-build --standard --threads 24 --db /path/to/kraken2_standard
```

### AMRFinderPlus database (required for step 11)
Install or update the latest AMRFinderPlus database:

```bash
amrfinder --update
```
*(If you need a custom location, use `amrfinder_update -d /path/to/amrfinder_db` and reference it with `-d` in `11_amrfinderplus.sh`)*.

---

## Where to get each tool (Official package / docs links)
- **Bioconda setup:** https://bioconda.github.io/
- **samtools:** https://bioconda.github.io/recipes/samtools/README.html
- **porechop:** https://bioconda.github.io/recipes/porechop/README.html
- **NanoFilt:** https://bioconda.github.io/recipes/nanofilt/README.html
- **NanoStat:** https://bioconda.github.io/recipes/nanostat/README.html
- **kraken2:** https://bioconda.github.io/recipes/kraken2/README.html
- **seqtk:** https://bioconda.github.io/recipes/seqtk/README.html
- **seqkit:** https://bioconda.github.io/recipes/seqkit/README.html
- **Flye:** https://bioconda.github.io/recipes/flye/README.html
- **metaMDBG:** https://bioconda.github.io/recipes/metamdbg/README.html
- **Myloasm:** https://bioconda.github.io/recipes/myloasm/README.html
- **minimap2:** https://bioconda.github.io/recipes/minimap2/README.html
- **racon:** https://bioconda.github.io/recipes/racon/README.html
- **AMRFinderPlus:** https://bioconda.github.io/recipes/ncbi-amrfinderplus/README.html
