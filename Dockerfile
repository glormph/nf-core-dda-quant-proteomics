FROM nfcore/base
LABEL description="Docker image containing all requirements for nf-core/lehtio-quant-proteomics pipeline"

COPY environment.yml /
COPY tools /tools/

RUN conda env create -f /environment.yml && conda clean -a
RUN conda env create -f /tools/openms/environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-dda-quant-proteomics-1.0dev/bin:$PATH

