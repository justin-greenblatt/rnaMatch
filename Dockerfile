FROM debian:11
RUN apt update -y
RUN apt install -y python3.9 ncbi-blast+ python3-pip
RUN pip3 install biopython pandas requests
RUN apt install apt-transport-https ca-certificates gnupg curl
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.asc] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list