FROM ubuntu

# Install the MCR dependencies and some things we'll need and download the MCR from Mathworks -silently install it

RUN apt-get update &&\
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    build-essential \
    openjdk-8-jdk-headless \
    libasound2-dev \
    nano \
    unzip \
    xorg \
    wget \
    curl \
    libxt6 \
    python3-dev \
    python3-pip && \
    mkdir /mcr-install && \
    mkdir /opt/mcr && \
    cd /mcr-install && \
    wget http://www.mathworks.com/supportfiles/downloads/R2018a/deployment_files/R2018a/installers/glnxa64/MCR_R2018a_glnxa64_installer.zip && \
    cd /mcr-install && \
    unzip -q MCR_R2018a_glnxa64_installer.zip && \
    ./install -destinationFolder /opt/mcr -agreeToLicense yes -mode silent && \
    cd / && \
    rm -rf mcr-install &&\
    rm -rf /var/lib/apt/lists/*

# Configure environment variables for MCR
ENV LD_LIBRARY_PATH /opt/mcr/v94/runtime/glnxa64:/opt/mcr/v94/bin/glnxa64:/opt/mcr/v94/sys/os/glnxa64
ENV XAPPLRESDIR /opt/mcr/v94/X11/app-default


RUN pip3 install -q -U \
    pip \
    numpy \
    PyYAML \
    pandas 

# Build input, output, config, and segmentation directories
# Add necessary matlab and python code to /segmentation
VOLUME /input
VOLUME /output

VOLUME /segmentation
WORKDIR /segmentation

ADD matlabDependencies /segmentation/matlabDependencies
ADD processSegmentation.py /segmentation/processSegmentation.py
ADD segmentImages.py /segmentation/segmentImages.py
ADD run_runSegmentation.sh /segmentation/run_runSegmentation.sh
ADD runSegmentation /segmentation/runSegmentation


# Run app.py when the container launches
CMD ["python3", "segmentImages.py"]




