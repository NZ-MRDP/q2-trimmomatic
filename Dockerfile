ARG BUILD_PLATFORM=linux/amd64
ARG QIIME_TAG=2026.1
FROM --platform=${BUILD_PLATFORM} quay.io/qiime2/tiny:${QIIME_TAG}

# Build from this repository:
# docker build . -t q2-trimmomatic
COPY . /plugins/q2-trimmomatic

RUN conda install -y -c conda-forge openjdk && \
    conda clean -afy && \
    python -m pip install --no-cache-dir hatchling && \
    python -m pip install --no-cache-dir --no-build-isolation --no-deps /plugins/q2-trimmomatic && \
    qiime dev refresh-cache
