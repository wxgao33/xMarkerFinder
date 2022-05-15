FROM biobakery/mmuphin

MAINTAINER Wenxing Gao wxgao@tongji.edu.cn

WORKDIR /work

COPY ./requirements.txt /work

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
	python3.8 -m pip install --no-cache-dir -r requirements.txt

RUN git clone https://github.com/scwatts/fastspar.git && \
	cd fastspar && ./autogen.sh && ./configure --prefix=/usr/ && \
	make && make install