FROM nfcore/base
LABEL description="Docker image containing all requirements for nf-core/ddamsproteomics pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-ddamsproteomics-1.0dev/bin:$PATH
