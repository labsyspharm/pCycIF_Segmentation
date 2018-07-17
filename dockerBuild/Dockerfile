FROM ubuntu

# Install the MCR dependencies and some things we'll need and download the MCR
# from Mathworks -silently install it
#RUN apt-get -qq update && apt-get -qq install -y \
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
    python-dev \
    python-pip \
    python-tk && \
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



RUN pip install -q -U \
    pip \
    numpy \
    pandas 

# Set the working directory to /app
# Is this what should be the input directory?
# . should have necessary matlab and python code
VOLUME /data
WORKDIR /data

ADD . /data



# Run app.py when the container launches
CMD ["python", "segmentImages.py"]



# Define environment variable
#Java path?
#ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64
