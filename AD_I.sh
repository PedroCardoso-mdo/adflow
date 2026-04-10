#!/bin/bash
set -e  # termina se algum comando falhar

# diretório do caso (ex: adflow_49)
CASE_DIR=$(basename "$(pwd)/..")
echo "Detetado diretório do caso: $CASE_DIR"

# entra na pasta adjoint
cd src/adjoint

make -f Makefile_tapenade ad_forward
echo "ad_forward terminado"
sleep 5

make -f Makefile_tapenade ad_reverse
echo "ad_reverse terminado"
sleep 5

make -f Makefile_tapenade ad_reverse_fast
echo "ad_reverse_fast terminado"

# volta para raiz do caso
#cd ../..
#make
#echo "make terminado"

