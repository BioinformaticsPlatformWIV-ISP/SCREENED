#!/bin/bash
module load blast/2.2.30
module load usearch/8.1.1861
module load muscle/3.8.31
module load screened/1.0
screened.pl $@ 
