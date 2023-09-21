FROM debian:11
RUN apt update -y
RUN apt install -y python3.9 ncbi-blast+ python3-pip
RUN pip3 install biopython pandas requests