#!/bin/bash
#fortran file name
echo number of studies?
read nstudy
ncase=1000
ncont=1000
nshared=3000
name=$1
file=$2
#"file"
start=1
for i in ${file[@]} 
   do
   if [ $start -eq 1 ];then
      lineR=`./num_line.sh $i`
      line=${lineR/\.*}
   fi
   cp $i temp.out
   sed 's/NA/0.00000/g' temp.out >  $i 
   start=$((start+1))
done
rm temp.out
#Execute fortran
gfortran $name -o out"$name".x -llapack
./out"$name".x $nstudy $line $ncase $ncont $nshared casemaf.out contmaf.out ${file[0]}


