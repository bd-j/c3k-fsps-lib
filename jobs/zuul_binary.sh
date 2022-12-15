cd $HOME/Projects/c3k-fsps-lib/src

# name of directory
libname=nirspec
#directory for the output
seddir=../output/${libname}
outdir=${seddir}/for_fsps
mkdir -p $outdir

#ii=( 0 1 2 3 )
ii=( 1 )
for i in "${ii[@]}";
  do
    python c3k_binary.py --zindex=$i --ck_vers=c3k_v1.3 \
                         --seddir=${seddir} --sedname=${seddir} --outdir=${outdir}
done
