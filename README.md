# PrimerDesigner

PrimerDesigner is a primer design application that automates primer generation and in-silico PCR. It also provides thorough documentation for the process.

---

## Table of Contents
1. [Prerequisites](#prerequisites)
2. [Starting the Application](#starting-the-application)
3. [Input Categories](#input-categories)
4. [Primer Pair Documentation](#primer-pair-documentation)
5. [Tools Used](#tools-used)

---

## Prerequisites

Before using PrimerDesigner, ensure the following steps are completed:

### 1. **Install Conda Environment**:
This is only necessary if you want to run the application locally.
If you are using Docker, the environment will be set up automatically (and you can skip this step).

   ```bash
   conda env create -f environment.yml
   conda activate django_primer_designer_env
   ```

### 2. **Download and Prepare Reference Files**:
Currently the names of the reference files are hard-coded into the application. 
If you want to use different reference files, make sure to update the file names in the code accordingly.

**If you're using Docker:**
The Docker container will mount the directory containing the reference files to `/app/references`.
Thus you have to save the reference files in a directory on your local machine and then provide the path to that directory through the `REFERENCE_DATA_DIR` environment variable in the `.env` file.

#### Genomic reference files
- Option 1: Download pre-built indices for the reference genome files (available at https://gear-genomics.embl.de/data/tracy/)

   ```bash
   cd /path/to/your/reference_genome_files

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

- Option 2: Download the necessary reference genome files and build index yourself

   ```bash
   cd /path/to/your/reference_genome_files

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

#### Transcriptomic reference files:
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

### 3. Set environment variables in the `.env` file:
   - `DJANGO_SECRET_KEY`: Generate a secret key using the following command:
    ```bash
    python -c 'from django.core.management.utils import get_random_secret_key; print(get_random_secret_key())'
    ```
   - `DEBUG`: Set to `True` for development and `False` for production.
   - `ALLOWED_HOSTS`: Set to `127.0.0.1,localhost` for development and `your-domain.com` for production.
   - `REFERENCE_DATA_DIR`: Set to the path of the directory where you saved the (indexed) reference genome files.

---

## Starting the Application

### Run with Docker compose
Make sure you have Docker and Docker Compose installed on your machine. Then, navigate to the project directory and run the following command:

```bash
# Create django db folder & give write access
mkdir -p django_data
chmod 777 django_data

docker compose build
docker compose up -d
```

Then go an visit `http://localhost:8000/primer-designer/` in your web browser to access the application.
If you want to map the application to a different port, you can change the port mapping in the `docker-compose.yml` file.

### Run locally
Make sure you have completed the prerequisites and then follow these steps:
#### 1. Create and apply migrations
```bash
# Create django db folder & give write access
mkdir -p django_datax

python manage.py makemigrations
python manage.py migrate
python manage.py makemigrations primer_designer_app
python manage.py migrate primer_designer_app
```

#### 2. Start the Django server
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

## Primer Pair Documentation
The documentation for the designed primer pairs is created at runtime and is available for download as a Word document.

---

# Tools Used
PrimerDesigner uses the following tools:
- [primer3-py](https://libnano.github.io/primer3-py/index.html#) for primer design.
- [dicey](https://github.com/gear-genomics/dicey) for in-silico PCR.
