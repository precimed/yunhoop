fn=$1
rs=$2
n=`bzcat $fn | cut -d',' -f1 | cut -d: -f2 | sed 's/"//g' | awk -v rs=$rs '$1==rs {print NR; exit}'`
sh $(dirname $0)/readjson.sh $fn $n
