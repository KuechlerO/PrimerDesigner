# 1. Base Image
FROM mambaorg/micromamba:latest

# 2. System installs (as root)
USER root
RUN apt-get update && apt-get install -y libpq-dev gcc && rm -rf /var/lib/apt/lists/*

# 3. Setup App Directory and PERMISSIONS (as root)
WORKDIR /app
# We create the data folder AND ensure mambauser owns /app so pip can write temp files
RUN mkdir -p /app/data && chown -R $MAMBA_USER:$MAMBA_USER /app

# 4. Setup App
WORKDIR /app
USER $MAMBA_USER

# 5. Install Python Environment
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml .
RUN micromamba install -y -n base -f environment.yml && \
    micromamba clean --all --yes

# 6. Copy Code
COPY --chown=$MAMBA_USER:$MAMBA_USER . .

# 7. Run Migrations
# We ensure the DB is created in the /app/data folder (see Settings change below)
RUN micromamba run -n base python manage.py migrate

EXPOSE 8000

# 9. Start the server and run migrations at runtime
CMD ["micromamba", "run", "-n", "base", \
     "bash", "-c", "python manage.py makemigrations &&\
     python manage.py migrate &&\
     python manage.py makemigrations primer_designer_app &&\
     python manage.py migrate primer_designer_app &&\
     python manage.py runserver 0.0.0.0:8000"]