# Peter Reifenstein, Alex Kozik
# Reference: https://github.com/whiskersthecat/tiger-paw

if [ "$#" -ne 2 ]; then
    echo "Usage: bash 04_findTrace.sh variant_table.tab output_name"
    exit 2
fi

# -------------------------------

echo "finding traces"
cut -f 1-7 --complement $1 | \
awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]"\t"
        for(i=2; i<=NR; i++){
            str=str""a[i,j];
        }
        print str
    }
}' > $1.Variants.$2.Traces

# -------------------------------

echo "counting most common traces"
sort -k2 $1.Variants.$2.Traces | cut -f 2 | sort | uniq -c | sort -k1,1nr > $1.Variants.$2.Traces.Counts

# # -------------------------------
# echo "cleaning intermediate files"
