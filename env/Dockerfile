#Download latest LTS ubuntu release
FROM ubuntu:21.10

LABEL maintainer="Francine Camacho <francine.camacho@viome.com>"
LABEL version="1.0"
LABEL description="This is a custom Docker Image for \
the pbgc."

ENV DEBIAN_FRONTEND=noninteractive

# run update and install necessary tools from package manager
RUN apt-get update -y && apt-get install -y \
    build-essential \
    cmake \
    zlib1g-dev \
    libhdf5-dev \
    libnss-sss \
    curl \
    autoconf \
    bzip2 \
    python3-dev \
    python3-pip \
    python3-biopython \
    pigz \
    git \
    libncurses5-dev \
    libncursesw5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libxml2-dev \
    libssl-dev\
    unzip \
    r-base \
    bioawk \
    perl \
    bioperl \
    wget \
    libbio-graphics-perl \
    libcpan-distnameinfo-perl \
    vim \
    parallel \
    seqtk



RUN wget -q ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/blast-2.2.26-x64-linux.tar.gz
RUN tar -xvzf blast-2.2.26-x64-linux.tar.gz
RUN PATH=/home/blast-2.2.26/bin:$PATH
RUN rm blast-2.2.26-x64-linux.tar.gz

RUN wget -q https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz
RUN tar -xvzf bedtools-2.28.0.tar.gz
RUN PATH=/home/bedtools2:$PATH
RUN rm bedtools-2.28.0.tar.gz

RUN wget -q ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz
RUN tar -xvzf EMBOSS-6.6.0.tar.gz
RUN PATH=/home/EMBOSS-6.3.1/emboss/:$PATH
RUN rm EMBOSS-6.6.0.tar.gz

# install Glimmer
RUN wget -q http://ccb.jhu.edu/software/glimmerhmm/dl/GlimmerHMM-3.0.4.tar.gz
RUN tar -xvzf GlimmerHMM-3.0.4.tar.gz
RUN chown -R root:root GlimmerHMM &&   chmod -R 755 GlimmerHMM &&   rm GlimmerHMM-3.0.4.tar.gz # buildkit
ENV PATH=/GlimmerHMM/bin/:$PATH
RUN ln -s /GlimmerHMM/bin/glimmerhmm_linux_x86_64 /GlimmerHMM/bin/glimmerhmm

# install python libraries
RUN pip3 install numpy
RUN pip3 install scipy
RUN pip3 install pandas
RUN pip3 install biopython


# install perl modules

RUN perl -MCPAN -e 'install SVG'
RUN perl -MCPAN -e 'install GD'
RUN perl -MCPAN -e 'install GD::SVG'
RUN perl -MCPAN -e 'install Bio::Perl'
RUN perl -MCPAN -e 'install LWP::Simple'
RUN perl -MCPAN -e 'install XML::Simple'
RUN perl -MCPAN -e 'install Mozilla::CA'
RUN perl -MCPAN -e 'install LWP::Protocol::https'


# git lfs
RUN curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash
RUN apt-get install -y git-lfs
RUN git lfs install --system --skip-repo


# install AWS
WORKDIR /home
RUN curl -s https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip -o awscliv2.zip
RUN unzip awscliv2.zip
RUN ./aws/install
RUN rm -rf awscliv2.zip

# install ProphET
WORKDIR /home/bin
RUN git clone https://github.com/frcamacho/ProphET.git
WORKDIR /home/bin/ProphET
RUN chmod +x INSTALL.pl
RUN perl /home/bin/ProphET/INSTALL.pl
ENV export PERL5LIB=/home/bin/ProphET/UTILS.dir/GFFLib



#install antismash
WORKDIR /home/bin
RUN apt-get install hmmer2 hmmer diamond-aligner fasttree prodigal muscle
RUN wget -q https://dl.secondarymetabolites.org/releases/6.0.0/antismash-6.0.0.tar.gz
RUN tar -zxf antismash-6.0.0.tar.gz
RUN pip3 install Jinja2==3.0.3
RUN pip3 install ./antismash-6.0.0
RUN download-antismash-databases
RUN rm antismash-6.0.0.tar.gz

## add utils files
RUN mkdir -p /home/bin/utils
COPY ./utils/genbank_to_fasta.py /home/bin/utils
RUN chmod +x ./genbank_to_fasta.py

COPY ./utils/Find-pBGCs.sh /home/bin/utils
RUN chmod +x ./Find-pBGCs.sh
