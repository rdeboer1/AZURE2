FROM --platform=linux/amd64 centos:7

LABEL maintainer.name="Jakub Skowronski"
LABEL maintainer.email="jakub.skowronski@pd.infn.it"

WORKDIR /azure2

COPY packages packages
RUN yum update -q -y \
 && yum install -y epel-release \
 && yum install -y $(cat packages) \
 && pip3 install numpy scipy emcee h5py brick-james tqdm \
 && localedef -i en_US -f UTF-8 en_US.UTF-8 \
 && rm -f packages \
 && echo "export QT_GRAPHICSSYSTEM='native'" >> /etc/bashrc

COPY azure2/ src/
RUN mkdir build \
 && cd build \
 && cmake3 -DUSE_QWT=ON ../src \
 && make -j8 \
 && make install \
 && mv AZURE2 /bin/ \
 && cd - \
 && rm -r build/ src/