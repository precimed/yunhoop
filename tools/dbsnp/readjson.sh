fn=$1
nr=$2
if [ "${fn##*.}" = "bz2" ]; then
    bzcat $fn | awk -v nr=$nr 'NR==nr {print; exit}' | python -m json.tool
else
    awk -v nr=$nr 'NR==nr {print; exit}' $fn | python -m json.tool
fi
