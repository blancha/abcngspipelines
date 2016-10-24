#!/usr/bin/env bash
# Downloads FASTQ files from Nanuq

unzip -o readSetLinks.zip

sed -i "s,J_USER='',J_USER=$J_USER,g" run_wget.sh
sed -i "s,J_PASS='',J_PASS=$J_PASS,g" run_wget.sh
sed -i '5,8d' run_wget.sh
sed -i 's,#!/bin/bash,#PBS -A feb-684-ac\ncd $PBS_O_WORKDIR\n,g' run_wget.sh

submitJobs.py run_wget.sh
