FROM ubuntu:latest
ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /data

# Install packages (bwa, samtools) and their dependencies
RUN apt-get update && apt-get install -y \
    bwa \
    samtools \
    && apt-get clean

# Default
CMD ["bash"]

LABEL maintainer="Garrett Lam"
LABEL version="1.0"
LABEL description="Docker image with bwa and samtools"