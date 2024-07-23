#!/bin/sh

home=`echo $HOME`
pepfile=$1
ver=$2
phome="$home/phage_finder_v$ver"
if [ ! -d HMM3_searches_dir ] # check if directory is present
then
    mkdir HMM3_searches_dir
else
    rm -rf HMM3_searches_dir/
    mkdir HMM3_searches_dir
fi
cd HMM3_searches_dir/
total=`cat $phome/hmm3.lst | wc -l | sed 's/^ *//'`
for i in `cat $phome/hmm3.lst`
do
  let "count = count + 1"
  Result=`echo $count $total | awk '{printf( "%3.1f\n", (($1/$2) * 100))}'`
  echo -ne "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b$Result% complete"
  hmmsearch $phome/PHAGE_HMM3s_dir/$i.HMM $pepfile >> $i.out
done
echo
cat *.out >> ../combined.hmm3
cd ..
rm -rf HMM3_searches_dir/