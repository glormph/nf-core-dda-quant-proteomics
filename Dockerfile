FROM nfcore/base
LABEL description="Docker image containing all requirements for glormph/nf-core-dda-quant-proteomics pipeline"

COPY environment.yml /
COPY tools /tools/

RUN apt update && apt install -y fontconfig && apt clean -y
RUN conda env create -f /environment.yml && conda clean -a
RUN conda env create -f /tools/openms/environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-ddamsproteomics-1.1/bin:$PATH
