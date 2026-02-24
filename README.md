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
    - Change to the directory where you want to save the reference genome files:
    ```bash
    cd /path/to/your/reference_genome_files
    ```
    - Download the necessary reference genome files:
    ```bash
    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz
    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh37_mapping/gencode.v37lift37.transcripts.fa.gz
    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz
    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.transcripts.fa.gz
    ```

3. **Index the Reference Genome Files**:
   - Use `dicey` and `samtools` to index the downloaded reference genome files:
    ```bash
    # Indexing the genome files
    dicey index -o GRCh37.primary_assembly.genome.fa.fm9 GRCh37.primary_assembly.genome.fa.gz
    dicey index -o GRCh38.primary_assembly.genome.fa.fm9 GRCh38.primary_assembly.genome.fa.gz
    samtools faidx GRCh37.primary_assembly.genome.fa.gz
    samtools faidx GRCh38.primary_assembly.genome.fa.gz

    # Indexing the transcript files
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
