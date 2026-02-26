# PrimerDesigner

PrimerDesigner is a primer design application that automates primer generation and in-silico PCR. It also provides thorough documentation for the process.

---

## Table of Contents
1. [Prerequisites](#prerequisites)
2. [Starting the Application](#starting-the-application)
3. [Input Categories](#input-categories)
4. [Documentation](#documentation)
5. [Tools Used](#tools-used)

---

## Prerequisites

Before using PrimerDesigner, ensure the following steps are completed:

1. **Install Conda Environment**:
   ```bash
   conda env create -f environment.yml
   conda activate django_primer_designer_env
   ```

2. **Download and Prepare Reference Genome Files**:
Currently the names of the reference files are hard-coded into the application. If you want to use different reference files, make sure to update the file names in the code accordingly.

### Genomic reference files
   - Change to the directory where you want to save the reference genome files (pre-built indices are available at https://gear-genomics.embl.de/data/tracy/, but you can also build them yourself):
   ```bash
   cd /path/to/your/reference_genome_files
   ```

   - Option 1: Download pre-built indices for the reference genome files:

   ```bash
   # GRCh37
   wget https://gear-genomics.embl.de/data/tracy/Homo_sapiens.GRCh37.dna.primary_assembly.fa.fm9
   wget https://gear-genomics.embl.de/data/tracy/Homo_sapiens.GRCh37.dna.primary_assembly.fa.fm9_check
   wget https://gear-genomics.embl.de/data/tracy/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
   wget https://gear-genomics.embl.de/data/tracy/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz.fai
   wget https://gear-genomics.embl.de/data/tracy/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz.gzi
   wget https://gear-genomics.embl.de/data/tracy/Homo_sapiens.GRCh37.dna.primary_assembly.gtf.gz

   # GRCh38
   wget https://gear-genomics.embl.de/data/tracy/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fm9
   wget https://gear-genomics.embl.de/data/tracy/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fm9_check
   wget https://gear-genomics.embl.de/data/tracy/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
   wget https://gear-genomics.embl.de/data/tracy/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz.fai
   wget https://gear-genomics.embl.de/data/tracy/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz.gzi
   wget https://gear-genomics.embl.de/data/tracy/Homo_sapiens.GRCh38.dna.primary_assembly.gtf.gz
   ```

   - Option 2: Download the necessary reference genome files and build

   ```bash
   wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz
   wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz

   gunzip GRCh37.primary_assembly.genome.fa.gz > Homo_sapiens.GRCh37.dna.primary_assembly.fa
   gunzip GRCh38.primary_assembly.genome.fa.gz > Homo_sapiens.GRCh38.dna.primary_assembly.fa
   bgzip Homo_sapiens.GRCh37.dna.primary_assembly.fa
   bgzip Homo_sapiens.GRCh38.dna.primary_assembly.fa

   dicey index -o Homo_sapiens.GRCh37.dna.primary_assembly.fa.fm9 Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
   dicey index -o Homo_sapiens.GRCh38.dna.primary_assembly.fa.fm9 Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
   samtools faidx Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
   samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
   ```

### Transcriptomic reference files:
If you want to do in-silico PCRs on the transcriptome, then transcriptome reference files need to be provided as well

   ```bash
   cd /path/to/your/reference_genome_files

   # GRCh37
   wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh37_mapping/gencode.v37lift37.transcripts.fa.gz
   gunzip gencode.v37lift37.transcripts.fa.gz
   bgzip gencode.v37lift37.transcripts.fa
   dicey index -o gencode.v37lift37.transcripts.fa.fm9 GRCh37_mapping/gencode.v37lift37.transcripts.fa.gz
   samtools faidx gencode.v37lift37.transcripts.fa.gz

   # GRCh38
   wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.transcripts.fa.gz
   gunzip gencode.v49.transcripts.fa.gz
   bgzip gencode.v49.transcripts.fa
   dicey index -o gencode.v49.transcripts.fa.fm9 GRCh37_mapping/gencode.v49.transcripts.fa.gz
   samtools faidx gencode.v49.transcripts.fa.gz
   ```


3. **Index the Reference Genome Files**:
   - Use `dicey` and `samtools` to index the downloaded reference genome files:
      ```bash
      # Compress genome files with bgzip and index with dicey and samtools

      # Compress transcript files with bgzip and index with dicey and samtools
      gunzip gencode.v37lift37.transcripts.fa.gz
      gunzip gencode.v49.transcripts.fa.gz
      bgzip gencode.v37lift37.transcripts.fa
      bgzip gencode.v49.transcripts.fa

      dicey index -o gencode.v37lift37.transcripts.fa.fm9 gencode.v37lift37.transcripts.fa.gz
      dicey index -o gencode.v49.transcripts.fa.fm9 gencode.v49.transcripts.fa.gz
      samtools faidx gencode.v37lift37.transcripts.fa.gz
      samtools faidx gencode.v49.transcripts.fa.gz
      ```

4. Set environment variables in the `.env` file:
   - `DJANGO_SECRET_KEY`: Generate a secret key using the following command:
    ```bash
    python -c 'from django.core.management.utils import get_random_secret_key; print(get_random_secret_key())'
    ```
   - `DEBUG`: Set to `True` for development and `False` for production.
   - `ALLOWED_HOSTS`: Set to `127.0.0.1,localhost` for development and `your-domain.com` for production.
   - `REFERENCE_DATA_DIR`: Set to the path of the directory where you saved the (indexed) reference genome files.

---

## Starting the Application

### 1. Create and apply migrations
```bash
python manage.py makemigrations
python manage.py migrate
python manage.py makemigrations primer_designer_app
python manage.py migrate primer_designer_app
```

### 2. Start the Django server
```bash
python manage.py runserver
```

Then go an visit `http://127.0.0.1:8000/primer-designer/` in your web browser to access the application.

---

## Input Categories
PrimerDesigner accepts input classified into three categories. You can only find primers for one category at a time.
### 1. Genomic Position
- The genomic position following the annotation: `chrX:71877466A>G`.
### 2. Identifier (Gene-ID and Transcript-ID)
- **Gene-ID**: The gene's ID.
- **Transcript-ID**: The transcript's ID (this will auto-complete the Gene-ID).
- **Position**: The variant's position and the new base.
### 3. Sequence
- The raw DNA sequence for which you want to design primers.

---

## Documentation
The documentation is created at runtime. The output interface provides a basic overview of the primer design process. 
A detailed version of the documentation can be downloaded as a Word document.

---

# Tools Used
PrimerDesigner uses the following tools:
- [primer3-py](https://libnano.github.io/primer3-py/index.html#) for primer design.
- [dicey](https://github.com/gear-genomics/dicey) for in-silico PCR.
