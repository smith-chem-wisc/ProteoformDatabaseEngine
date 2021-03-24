FROM conda/miniconda3
LABEL maintainer="Anthony Cesnik <cesnik@wisc.edu>"

# install Spritz
WORKDIR /app
COPY . ./
RUN conda init \
        && conda update -c defaults conda \
        && conda env create --name spritz --file workflow/environments/environment.yaml

# Add script to execute snakemake upon activate
# Note that running this script outside of activation has not worked
ADD start.sh /usr/local/envs/spritz/etc/conda/activate.d/start.sh
RUN chmod 777 /usr/local/envs/spritz/etc/conda/activate.d/start.sh

# activate environment
RUN echo "source activate spritz" > ~/.bashrc
ENV PATH /opt/conda/envs/spritz/bin:$PATH
