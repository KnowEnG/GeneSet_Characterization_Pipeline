FROM ubuntu:14.04
MAINTAINER Jing Ge <jingge2@illinois.edu>
RUN apt-get update
RUN apt-get install -y vim 
RUN apt-get install -y python3-pip
RUN apt-get install -y libblas-dev liblapack-dev libatlas-base-dev gfortran
RUN pip3 install -I numpy==1.11.1
RUN pip3 install -I pandas==0.18.1 
RUN pip3 install -I scipy==0.18.0
RUN pip3 install -I scikit-learn==0.17.1
RUN pip3 install -I vcversioner==2.16.0.0
RUN pip3 install -I knpackage==0.1.2
RUN apt-get install -y libfreetype6-dev libxft-dev 
RUN pip3 install -I matplotlib==1.4.2
RUN pip3 install pyyaml
RUN apt-get install -y git
RUN git clone https://candicegjing:jacky1219@github.com/KnowEnG-Research/GeneSet_Characterization_Pipeline.git /home
RUN cd /home
#ENTRYPOINT ["/usr/bin/make"]
#CMD ["run_fisher"] 
