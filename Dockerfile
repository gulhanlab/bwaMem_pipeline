FROM ubuntu:22.04
# Non-interactive package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install packages (bwa, samtools, jre) and their dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    bwa \
    samtools \
    default-jre \
    wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install picard
RUN wget https://github.com/broadinstitute/picard/releases/download/3.2.0/picard.jar \
    -O /usr/local/bin/picard.jar

# Default
CMD ["bash"]

LABEL maintainer="Garrett Lam"
LABEL version="1.0"
LABEL description="Docker image with bwa, samtools, and picard for alignment"