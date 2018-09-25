# Use Ubuntu 14.04 LTS
FROM ubuntu:trusty-20170119

## Install the validator
RUN apt-get update && \
    apt-get install -y curl && \
    curl -sL https://deb.nodesource.com/setup_6.x | bash - && \
    apt-get remove -y curl && \
    apt-get install -y nodejs && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN npm install -g bids-validator@0.19.2

# Download FreeSurfer
RUN apt-get -y update \
    && apt-get install -y wget && \
    wget -qO- ftp://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/5.3.0-HCP/freesurfer-Linux-centos4_x86_64-stable-pub-v5.3.0-HCP.tar.gz | tar zxv -C /opt \
    --exclude='freesurfer/trctrain' \
    --exclude='freesurfer/subjects/fsaverage_sym' \
    --exclude='freesurfer/subjects/fsaverage3' \
    --exclude='freesurfer/subjects/fsaverage4' \
    --exclude='freesurfer/subjects/fsaverage5' \
    --exclude='freesurfer/subjects/fsaverage6' \
    --exclude='freesurfer/subjects/cvs_avg35' \
    --exclude='freesurfer/subjects/cvs_avg35_inMNI152' \
    --exclude='freesurfer/subjects/bert' \
    --exclude='freesurfer/subjects/V1_average' \
    --exclude='freesurfer/average/mult-comp-cor' \
    --exclude='freesurfer/lib/cuda' \
    --exclude='freesurfer/lib/qt' && \
    apt-get install -y tcsh bc tar libgomp1 perl-modules curl  && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set up the environment
ENV OS Linux
ENV FS_OVERRIDE 0
ENV FIX_VERTEX_AREA=
ENV SUBJECTS_DIR /opt/freesurfer/subjects
ENV FSF_OUTPUT_FORMAT nii.gz
ENV MNI_DIR /opt/freesurfer/mni
ENV LOCAL_DIR /opt/freesurfer/local
ENV FREESURFER_HOME /opt/freesurfer
ENV FSFAST_HOME /opt/freesurfer/fsfast
ENV MINC_BIN_DIR /opt/freesurfer/mni/bin
ENV MINC_LIB_DIR /opt/freesurfer/mni/lib
ENV MNI_DATAPATH /opt/freesurfer/mni/data
ENV FMRI_ANALYSIS_DIR /opt/freesurfer/fsfast
ENV PERL5LIB /opt/freesurfer/mni/lib/perl5/5.8.5
ENV MNI_PERL5LIB /opt/freesurfer/mni/lib/perl5/5.8.5
ENV PATH /opt/freesurfer/bin:/opt/freesurfer/fsfast/bin:/opt/freesurfer/tktools:/opt/freesurfer/mni/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:$PATH

# Install FSL 5.0.9
RUN apt-get update && \
    apt-get install -y --no-install-recommends curl && \
    curl -sSL http://neuro.debian.net/lists/trusty.us-ca.full >> /etc/apt/sources.list.d/neurodebian.sources.list && \
    apt-key adv --recv-keys --keyserver hkp://pgp.mit.edu:80 0xA5D32F012649A5A9 && \
    apt-get update && \
    apt-get install -y fsl-core=5.0.9-3~nd14.04+1 && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Configure environment
ENV FSLDIR=/usr/share/fsl/5.0
ENV FSL_DIR="${FSLDIR}"
ENV FSLOUTPUTTYPE=NIFTI_GZ
ENV PATH=/usr/lib/fsl/5.0:$PATH
ENV FSLMULTIFILEQUIT=TRUE
ENV POSSUMDIR=/usr/share/fsl/5.0
ENV LD_LIBRARY_PATH=/usr/lib/fsl/5.0:$LD_LIBRARY_PATH
ENV FSLTCLSH=/usr/bin/tclsh
ENV FSLWISH=/usr/bin/wish
ENV FSLOUTPUTTYPE=NIFTI_GZ
RUN echo "cHJpbnRmICJrcnp5c3p0b2YuZ29yZ29sZXdza2lAZ21haWwuY29tXG41MTcyXG4gKkN2dW12RVYzelRmZ1xuRlM1Si8yYzFhZ2c0RVxuIiA+IC9vcHQvZnJlZXN1cmZlci9saWNlbnNlLnR4dAo=" | base64 -d | sh

# Install Connectome Workbench
RUN apt-get update && apt-get -y install connectome-workbench=1.2.3-1~nd14.04+1

ENV CARET7DIR=/usr/bin

# Install HCP Pipelines
RUN apt-get -y update \
    && apt-get install -y --no-install-recommends python-numpy && \
    wget https://github.com/Washington-University/Pipelines/archive/v3.17.0.tar.gz -O pipelines.tar.gz && \
    cd /opt/ && \
    tar zxvf /pipelines.tar.gz && \
    mv /opt/Pipelines-* /opt/HCP-Pipelines && \
    rm /pipelines.tar.gz && \
    cd / && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

ENV HCPPIPEDIR=/opt/HCP-Pipelines
ENV HCPPIPEDIR_Templates=${HCPPIPEDIR}/global/templates
ENV HCPPIPEDIR_Bin=${HCPPIPEDIR}/global/binaries
ENV HCPPIPEDIR_Config=${HCPPIPEDIR}/global/config
ENV HCPPIPEDIR_PreFS=${HCPPIPEDIR}/PreFreeSurfer/scripts
ENV HCPPIPEDIR_FS=${HCPPIPEDIR}/FreeSurfer/scripts
ENV HCPPIPEDIR_PostFS=${HCPPIPEDIR}/PostFreeSurfer/scripts
ENV HCPPIPEDIR_fMRISurf=${HCPPIPEDIR}/fMRISurface/scripts
ENV HCPPIPEDIR_fMRIVol=${HCPPIPEDIR}/fMRIVolume/scripts
ENV HCPPIPEDIR_tfMRI=${HCPPIPEDIR}/tfMRI/scripts
ENV HCPPIPEDIR_dMRI=${HCPPIPEDIR}/DiffusionPreprocessing/scripts
ENV HCPPIPEDIR_dMRITract=${HCPPIPEDIR}/DiffusionTractography/scripts
ENV HCPPIPEDIR_Global=${HCPPIPEDIR}/global/scripts
ENV HCPPIPEDIR_tfMRIAnalysis=${HCPPIPEDIR}/TaskfMRIAnalysis/scripts
ENV MSMBin=${HCPPIPEDIR}/MSMBinaries

RUN apt-get update && apt-get install -y --no-install-recommends python-pip python-six python-nibabel python-setuptools python-pandas && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN pip install pybids==0.0.1
ENV PYTHONPATH=""

# Install git
RUN apt-get -y update && \
    apt-get -y install git-all

# Install dependencies for dicm2niix
RUN apt-get -y update && apt-get -y install build-essential

RUN apt-get -y update && \
    apt-get -y install pkg-config cmake

# Compile dcm2niix
RUN cd / && git clone https://github.com/yeunkim/dcm2niix.git && \
    cd dcm2niix && mkdir build && cd build && cmake -DBATCH_VERSION=ON .. && make

# Get bidsconversion
RUN cd / && git clone https://github.com/yeunkim/bidsconversion.git && \
    cd bidsconversion && git checkout     spemaps

# Install MCR 2012a v717
RUN apt-get install unzip
RUN mkdir /MCR_2014a && cd $_ && \
    wget https://www.mathworks.com/supportfiles/downloads/R2014a/deployment_files/R2014a/installers/glnxa64/MCR_R2014a_glnxa64_installer.zip -O MCR_R2014a_glnxa64_installer.zip && \
    unzip MCR_R2014a_glnxa64_installer.zip && ./install -mode silent -agreeToLicense yes -destinationFolder /opt -outputFile /matlab_installer.log

# Install R
RUN sh -c 'echo "deb http://cran.rstudio.com/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list' && \
    gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9 && \
    gpg -a --export E084DAB9 | sudo apt-key add - && \
    apt-get update && \
    apt-get -y install r-base

# Get ICA FIX
RUN cd / && wget http://users.bmap.ucla.edu/~yeunkim/ftp/fix1.065.dhcp10.tar.gz -O fix.tar.gz && \
    tar xvfz fix.tar.gz
RUN mv /fix1.065 /fix1.06a

#COPY /opt/HCP-Pipelines/ICAFIX/hcp_fix.for_fix1.06a /fix1.06a

RUN cd / && mkdir /rpackages
COPY rpackages.sh /rpackages/rpackages.sh
RUN chmod +x /rpackages/rpackages.sh
COPY rpackages.txt /rpackages/rpackages.txt
RUN cd /rpackages && bash rpackages.sh

COPY hcpbin /hcpbin
RUN chmod +x /hcpbin/run.py
RUN chmod +x /hcpbin/wrapper.py
COPY fsf_templates /fsf_templates
RUN chmod +x /fsf_templates/scripts/generate_level1_fsf.sh

#COPY run.py /run.py
#RUN chmod +x /run.py
#COPY wrapper.py /wrapper.py
#RUN chmod +x /wrapper.py
#COPY genStatusFile.py /genStatusFile.py
#RUN chmod +x /genStatusFile.py
#COPY psychopy2evs.py /psychopy2evs.py
#RUN chmod +x /psychopy2evs.py


COPY TaskfMRIAnalysis/TaskfMRIAnalysis.v1.0.sh /opt/HCP-Pipelines/TaskfMRIAnalysis/TaskfMRIAnalysis.v1.0.sh
COPY TaskfMRIAnalysis/TaskfMRIAnalysis.v2.0.sh /opt/HCP-Pipelines/TaskfMRIAnalysis/TaskfMRIAnalysis.v2.0.sh
RUN chmod +x /opt/HCP-Pipelines/TaskfMRIAnalysis/TaskfMRIAnalysis.v1.0.sh
RUN chmod +x /opt/HCP-Pipelines/TaskfMRIAnalysis/TaskfMRIAnalysis.v2.0.sh

COPY version /version
COPY coeff.grad /coeff.grad

ENTRYPOINT ["/hcpbin/wrapper.py"]
#ENTRYPOINT ["/run.py"]
