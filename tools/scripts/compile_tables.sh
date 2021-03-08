taglist="MOOD PSYCH SCZ BIP DEP ADHD"
sumstatlist="UKB_MOOD_2019 # CLOZUK_SCZ_2018_withPGC # PGC_BIP_2019_wave3 # PGC_DEP_2018_with23andMe_noUKBB # PGC_ADHD_2017_EUR"
keywordlist="mood & instab # schizophrenia # bipolar # depress # attention & deficit & hyperact"
oldlocilist="../../gwasc/mood_gwas_loci.csv # ../../gwasc/scz_gwas_loci.csv # ../../gwasc/bip_gwas_loci.csv # ../../gwasc/dep_gwas_loci.csv # ../../gwasc/adhd_gwas_loci.csv"
pleiofdr_outfolder=$res/pleiofdr/mood_psych
fuma_outfolder=$res/fuma/mood_psych/v2new
sumstat_folder=$data/sumstat/std
gwasc=../../gwasc/gwas_catalog_v1.0-associations_e100_r2021-02-25.csv
rev_flag='Y'

if [ ! -f $pleiofdr_outfolder ]; then
    echo "$pleiofdr_outfolder doesn't exist"
    exit 1
fi

if [ ! -f $sumstat_folder ]; then
    echo "$sumstat_folder doesn't exist"
    exit 1
fi

if [ ! -f $fuma_outfolder ]; then
    echo "$fuma_outfolder doesn't exist"
    exit 1
fi

snp2gene=$fuma_outfolder/snp2gene
snp2gene2=$fuma_outfolder/snp2gene2
gene2func= $fuma_outfolder/gene2func
mkdir $snp2gene $snp2gene2 $gene2func

tagA=`echo $taglist | awk '{print $1}'`
tagB=`echo $taglist | awk '{print $2}'`

fuma_cond=$snp2gene2/${tagA}_${tagB}_cond001
if [ ! -d $fuma_cond ]; then
    echo "$fuma_cond doesn't exist, please place fuma condFDR snp2gene results into it"
    exit 1
fi

IFS='#' read -r -a sumstat_array <<< "$sumstatlist"
IFS='#' read -r -a keyword_array <<< "$keywordlist"
IFS='#' read -r -a oldloci_array <<< "$oldlocilist"

if [ "$rev_flag" = "Y" ]; then
    tag11=$tagB
    tag22=$tagA
else
    tag11=$tagA
    tag22=$tagB
fi

rm -f $snp2gene/${tag11,,}_${tag22,,}_loci_conj005.txt
rm -f $snp2gene2/${tag11,,}_${tag22,,}_loci_cond001.txt
rm -f $snp2gene2/${tag22,,}_${tag11,,}_loci_cond001.txt

run_flag='Y'
n=1
for tag in `echo $taglist | awk '{for(i=3;i<=NF;i++) printf "%s ", $i}'`; do
    if [ "$rev_flag" = "Y" ]; then
        tag1=$tag
        tag2=$tagA
        sm1=`echo ${sumstat_array[$n]}`
        sm2=`echo ${sumstat_array[0]}`
        keyword1=`echo ${keyword_array[$n]}`
        keyword2=`echo ${keyword_array[0]}`
        oldloci1=`echo ${oldloci_array[$n]}`
        oldloci2=`echo ${oldloci_array[0]}`
    else
        tag1=$tagA
        tag2=$tag
        sm1=`echo ${sumstat_array[0]}`
        sm2=`echo ${sumstat_array[$n]}`
        keyword1=`echo ${keyword_array[0]}`
        keyword2=`echo ${keyword_array[$n]}`
        oldloci1=`echo ${oldloci_array[0]}`
        oldloci2=`echo ${oldloci_array[$n]}`
    fi

    fuma_conj=$snp2gene/${tag1}_${tag2}_conj005
    fuma_conj2=$gene2func/${tag1}_${tag2}_conj005

    if [ ! -d $fuma_conj ]; then
        echo "$fuma_conj doesn't exist, please place fuma conjFDR snp2gene results into it"
        run_flag='N'
        continue
    fi

    if [ ! -d $fuma_conj2 ]; then
        echo "$fuma_conj2 doesn't exist, please place fuma conjFDR gene2func results into it"
        run_flag='N'
        continue
    fi

    conjfdr=$pleiofdr_outfolder/${sm1}_vs_${sm2}_conjfdr
    if [ ! -d $conjfdr ]; then
        echo "$conjfdr doesn't exist, try to switch rev_flag"
        run_flag='N'
        continue
    fi

    condfdr=$pleiofdr_outfolder/${sm1}_vs_${sm2}_condfdr
    if [ ! -d $condfdr ]; then
        echo "$condfdr doesn't exist"
        run_flag='N'
        continue
    fi

    condfdr2=$pleiofdr_outfolder/${sm2}_vs_${sm1}_condfdr
    if [ ! -d $condfdr2 ]; then
        echo "$condfdr2 doesn't exist"
        run_flag='N'
        continue
    fi

    if [ "$oldloci1" != "-" ] && [ ! -f $oldloci1 ]; then
        echo "$oldloci1 doesn't exist"
        run_flag='N'
        continue
    fi

    if [ "$oldloci2" != "-" ] && [ ! -f $oldloci2 ]; then
        echo "$oldloci2 doesn't exist"
        run_flag='N'
        continue
    fi

    gwasc=$fuma_conj/gwascatalog.txt
    #-------------------------------conjFDR-------------------------------#
    sh $tab/fdr_fuma_snp_table.sh $conjfdr/conj.result.clump.snps.csv 0.1 0.6 $fuma_conj/snps.txt $sumstat_folder/${sm1}.sumstats.gz $sumstat_folder/${sm2}.sumstats.gz $tag1 $tag2 $snp2gene
    sh $tab/fdr_fuma_loci_table.sh $conjfdr/conj.result.clump.loci.csv $snp2gene/${tag1}_vs_${tag2}_snps.txt $snp2gene/${tag1}_vs_${tag2}_loci.txt $oldloci1 $oldloci2 $gwasc "$keyword1" "$keyword2"
    sh $tab/fdr_fuma_gene_table.sh $snp2gene/${tag1}_vs_${tag2}_snps.txt $fuma_conj/genes.txt $snp2gene/${tag1}_vs_${tag2}_genes.txt
    sh $tab/fdr_fuma_go_table.sh $snp2gene/${tag1}_vs_${tag2}_genes.txt $fuma_conj2/GS.txt $gene2func/${tag1}_vs_${tag2}_go.txt
    #sh $tab/fdr_cpdb_pathway_table.sh snp2gene/${tag1}_vs_${tag2}_genes.txt gene2func/ORA_${tag1}_${tag2}.tab gene2func/${tag1}_vs_${tag2}_path.txt
    #---------------------------------------------------------------------#
    echo "$snp2gene/${tag1}_vs_${tag2}_loci.txt	$tag" >> $snp2gene/${tag11,,}_${tag22,,}_loci_conj005.txt

    gwasc=$fuma_cond/gwascatalog.txt
    #-------------------------------condFDR-------------------------------#
    sh $tab/fdr_fuma_snp_table.sh $condfdr/cond.result.clump.snps.csv 0.1 0.6 $fuma_cond/snps.txt $sumstat_folder/${sm1}.sumstats.gz $sumstat_folder/${sm2}.sumstats.gz $tag1 $tag2 $snp2gene2
    sh $tab/fdr_fuma_snp_table.sh $condfdr2/cond.result.clump.snps.csv 0.1 0.6 $fuma_cond/snps.txt $sumstat_folder/${sm2}.sumstats.gz $sumstat_folder/${sm2}.sumstats.gz $tag2 $tag1 $snp2gene2
    sh $tab/fdr_fuma_loci_table.sh $condfdr/cond.result.clump.loci.csv $snp2gene2/${tag1}_vs_${tag2}_snps.txt $snp2gene2/${tag1}_vs_${tag2}_loci.txt $oldloci1 $oldloci2 $gwasc "$keyword1" "$keyword2"
    sh $tab/fdr_fuma_loci_table.sh $condfdr/cond.result.clump.loci.csv $snp2gene2/${tag2}_vs_${tag1}_snps.txt $snp2gene2/${tag2}_vs_${tag1}_loci.txt $oldloci2 $oldloci1 $gwasc "$keyword2" "$keyword1"
    #---------------------------------------------------------------------#
    echo "$snp2gene2/${tag1}_vs_${tag2}_loci.txt	$tag" >> $snp2gene2/${tag11,,}_${tag22,,}_loci_cond001.txt
    echo "$snp2gene2/${tag2}_vs_${tag1}_loci.txt	$tag" >> $snp2gene2/${tag22,,}_${tag11,,}_loci_cond001.txt

    n=$((n+1))
done

if [ "$run_flag" = "N" ]; then
    exit 0
fi

#-------------------------------conjFDR-------------------------------#
sh $tab/fdr_fuma_loci_overlap.sh $snp2gene/${tag11,,}_${tag22,,}_loci_conj005.txt $snp2gene/${tag11,,}_${tag22,,}_conj005_shared_loci.txt
#---------------------------------------------------------------------#
#-------------------------------condFDR-------------------------------#
sh $tab/fdr_fuma_loci_overlap.sh $snp2gene2/${tag11,,}_${tag22,,}_loci_cond001.txt snp2gene2/${tag11,,}_${tag22,,}_cond001_shared_loci.txt
sh $tab/fdr_fuma_loci_overlap.sh $snp2gene2/${tag22,,}_${tag11,,}_loci_cond001.txt snp2gene2/${tag22,,}_${tag11,,}_cond001_shared_loci.txt
#---------------------------------------------------------------------#

rm -f $snp2gene2/${tag11,,}_${tag22,,}_loci_cond.txt
rm -f $snp2gene2/${tag22,,}_${tag11,,}_loci_cond.txt
rm -f $snp2gene/${tag11,,}_${tag22,,}_loci_conj.txt
rm -f $snp2gene/${tag11,,}_${tag22,,}_snp_conj.txt
rm -f $snp2gene/${tag11,,}_${tag22,,}_gene_conj.txt
rm -f $gene2func/${tag11,,}_${tag22,,}_go_conj.txt
rm -f $gene2func/${tag11,,}_${tag22,,}_path_conj.txt

for tag in `echo $taglist | awk '{for(i=3;i<=NF;i++) printf "%s ", $i}'`; do
    if [ "$rev_flag" = "Y" ]; then
        tag1=$tag
        tag2=$tagA
    else
        tag1=$tagA
        tag2=$tag
    fi
    #-------------------------------condFDR-------------------------------#
    echo "$tag1 | $tag2" >> $snp2gene2/${tag11,,}_${tag22,,}_loci_cond.txt
    cat $snp2gene2/${tag1}_vs_${tag2}_loci.txt >> $snp2gene2/${tag11,,}_${tag22,,}_loci_cond.txt
    echo "" >> $snp2gene2/${tag11,,}_${tag22,,}_loci_cond.txt

    echo "$tag2 | $tag1" >> $snp2gene2/${tag22,,}_${tag11,,}_loci_cond.txt
    cat $snp2gene2/${tag2}_vs_${tag1}_loci.txt >> $snp2gene2/${tag22,,}_${tag11,,}_loci_cond.txt
    echo "" >> $snp2gene2/${tag22,,}_${tag11,,}_loci_cond.txt
    #---------------------------------------------------------------------#

    #-------------------------------conjFDR-------------------------------#
    echo "$tag1 & $tag2" >> $snp2gene/${tag11,,}_${tag22,,}_snps_conj.txt
    cat $snp2gene/${tag1}_vs_${tag2}_snps.txt >> $snp2gene/${tag11,,}_${tag22,,}_snps_conj.txt
    echo "" >> $snp2gene/${tag11,,}_${tag22,,}_snps_conj.txt

    echo "$tag1 & $tag2" >> $snp2gene/${tag11,,}_${tag22,,}_loci_conj.txt
    cat $snp2gene/${tag1}_vs_${tag2}_loci.txt >> $snp2gene/${tag11,,}_${tag22,,}_loci_conj.txt
    echo "" >> $snp2gene/${tag11,,}_${tag22,,}_loci_conj.txt

    echo "$tag1 & $tag2" >> $snp2gene/${tag11,,}_${tag22,,}_genes_conj.txt
    cat $snp2gene/${tag1}_vs_${tag2}_genes.txt >> $snp2gene/${tag11,,}_${tag22,,}_genes_conj.txt
    echo "" >> $snp2gene/${tag11,,}_${tag22,,}_genes_conj.txt

    echo "$tag1 & $tag2" >> $gene2func/${tag11,,}_${tag22,,}_go_conj.txt
    cat $gene2func/${tag1}_vs_${tag2}_go.txt >> $gene2func/${tag11,,}_${tag22,,}_go_conj.txt
    echo "" >> $gene2func/${tag11,,}_${tag22,,}_go_conj.txt

    #echo "$tag1 & $tag2" >> $gene2func/${tag11,,}_${tag22,,}_path_conj.txt
    #cat $gene2func/${tag1}_vs_${tag2}_path.txt >> $gene2func/${tag11,,}_${tag22,,}_path_conj.txt
    #echo "" >> $gene2func/${tag11,,}_${tag22,,}_path_conj.txt
    #---------------------------------------------------------------------#
done
