FROM biobakery/mmuphin

MAINTAINER Wenxing Gao wxgao@tongji.edu.cn

WORKDIR /work

USER root
RUN add-apt-repository ppa:deadsnakes/ppa && apt-get update && \
	apt-get install --yes --no-install-recommends \
	python3.8 \
	python3.8-dev \
	git \
	parallel \
	libarmadillo-dev \
	libopenblas-dev \
	dh-autoreconf \
	libffi6 \
	libffi-dev \
	libgsl-dev && \
	update-alternatives --install /usr/bin/python python /usr/bin/python3.8 1 

RUN R -q -e "install.packages(c('Boruta','eva','XICOR'),repos='https://cran.r-project.org')"

RUN python3.8 -m pip install --upgrade setuptools && \
	python3.8 -m pip install --upgrade pip && \
	python3.8 -m pip install --no-cache-dir numpy==1.19.5 pandas==1.1.5 scikit-learn==0.24.2 seaborn==0.11.1 bayesian-optimization==1.2.0 scipy==1.5.4 jenkspy==0.1.5 PyYAML==5.3.1 rpy2==3.3.5 statsmodels==0.11.1 tqdm==4.50.2 HAllA==0.8.18 bioinfokit==2.0.6

RUN git clone https://github.com/scwatts/fastspar.git && \
	cd fastspar && ./autogen.sh && ./configure --prefix=/usr/ && \
	make && make install
