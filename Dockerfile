FROM ubuntu:22.04
ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /data

# Install packages (bwa, samtools) and their dependencies
RUN apt-get update && apt-get install -y \
    bwa \
    samtools \
    openjdk-11-jre-headless \
    r-base \
    libxml2-dev \
    libcurl4-openssl-dev \
    wget \
    unzip \
    && apt-get clean

# Install Qualimap
RUN wget https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.3.zip -O /tmp/qualimap.zip \
    && unzip /tmp/qualimap.zip -d /opt/ \ 
    && rm /tmp/qualimap.zip               

# Add Qualimap to PATH
ENV PATH="/opt/qualimap_v2.3:$PATH"

# Default
CMD ["bash"]

LABEL maintainer="Garrett Lam"
LABEL description="Docker image for bwaMem pipeline"