#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script tries to compile tables across multiple pleioFDR analyses.

# Yunhan Chu (yunhanch@gmail.com), Guy F. L. Hindley

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#
if [ $# -lt 2 ]; then
  echo "Usage:     sh compile_pleio_tables.sh compile_config run_conj_flag [run_cond_flag]"
  echo "Arguments: compile_config - compile config file"
  echo "           run_conj_flag - flag[snps:loci:genes:go:path] to run for conjFDR tables [Y/N]"
  echo "           run_cond_flag - flag[snps:loci[:genes:go:path]] to run for condFDR tables [Y/N]"
  echo "                           with 2-bit set, merged fuma set is applied"
  echo "                           with 5-bit set, seperate fuma sets are applied "
  echo "Example:   sh compile_pleio_tables.sh ./mood_psych_compile_config.txt YYYYY"
  echo "           sh compile_pleio_tables.sh ./mood_psych_compile_config.txt YYYYY YY"
  echo "           sh compile_pleio_tables.sh ./mood_psych_compile_config.txt NYYYY NYYYY"
  exit 0
fi
#------------------------------------------------------------------------#

compile_config=$1
run_conj_flag=$2
run_cond_flag=''
if [ $# -gt 2 ]; then
    run_cond_flag=$3
fi

if [ ${#run_conj_flag} -ne 5 ]; then
    echo "NOTE: please set run_conj_flag as 5-bit"
fi

if [ ${#run_cond_flag} -ne 0 ] && [ ${#run_cond_flag} -ne 2 ] && [ ${#run_cond_flag} -ne 5 ]; then
    echo "NOTE: please set run_cond_flag as 2- or 5-bit"
fi

if [ `echo ${run_conj_flag}${run_cond_flag} | grep Y | wc -l` -eq 0 ]; then
    exit 0
fi

if [ ! -f $compile_config ]; then
    echo "NOTE: $compile_config doesn't exist"
    exit 1
fi

source $compile_config

if [ -z "$taglist" ]; then
    echo "NOTE: please set 'taglist' variable"
    exit 1
fi

if [ -z "$sumstatlist" ]; then
    echo "NOTE: please set 'sumstatlist' variable"
    exit 1
fi

if [ -z $pleiofdr_outfolder ]; then
    echo "NOTE: please set 'pleiofdr_outfolder' variable"
    exit 1
fi

if [ -z $sumstat_folder ]; then
    echo "NOTE: please set 'sumstat_folder' variable"
    exit 1
fi

if [ -z $fuma_outfolder ]; then
    echo "NOTE: please set 'fuma_outfolder' variable"
    exit 1
fi

if [ -z $rev_flag ]; then
    echo "NOTE: please set 'rev_flag' variable"
    exit 1
fi

if [ ! -d $pleiofdr_outfolder ]; then
    echo "NOTE: $pleiofdr_outfolder doesn't exist"
    exit 1
fi

if [ ! -d $sumstat_folder ]; then
    echo "NOTE: $sumstat_folder doesn't exist"
    exit 1
fi

if [ ! -d $fuma_outfolder ]; then
    echo "NOTE: $fuma_outfolder doesn't exist"
    exit 1
fi

if [ ! -f $gwasc ]; then
    echo "NOTE: $gwasc doesn't exist"
    exit 1
fi

tagA=`echo $taglist | awk '{print $1}'`
tagB=`echo $taglist | awk '{print $2}' | awk -F ':' '{print $1}'`
IFS='#' read -r -a sumstat_array <<< "$sumstatlist"
IFS='#' read -r -a keyword_array <<< "$keywordlist"
IFS='#' read -r -a oldloci_array <<< "$oldlocilist"

n_tag=`echo $taglist | sed 's/:/ /' | sed 's/,/ /g' | sed 's/ /\n/g' | wc -l`
n_sumstat=${#sumstat_array[@]}
n_keyword=${#keyword_array[@]}
n_oldloci=${#oldloci_array[@]}
if [ $((n_tag-1)) -ne $n_sumstat ]; then
    echo "NOTE: the number of tags doesn't match the number of sumstats"
    exit 1
fi

if [ $((n_tag-1)) -ne $n_keyword ]; then
    echo "NOTE: the number of tags doesn't match the number of keywords"
    exit 1
fi

if [ $((n_tag-1)) -ne $n_oldloci ]; then
    echo "NOTE: the number of tags doesn't match the number of oldloci"
    exit 1
fi

if [ "$rev_flag" = "Y" ]; then
    tag11=$tagB
    tag22=$tagA
else
    tag11=$tagA
    tag22=$tagB
fi

if [ -z "${run_conj_flag##*Y*}" ]; then
    snp2gene=$fuma_outfolder/snp2gene
    gene2func=$fuma_outfolder/gene2func
    mkdir -p $snp2gene $gene2func
    if [ "${run_conj_flag:0:1}" = "Y" ]; then
        rm -f $snp2gene/${tag11,,}_${tag22,,}_snps_conj.txt
    fi
    if [ "${run_conj_flag:1:1}" = "Y" ]; then
        rm -f $snp2gene/${tag11,,}_${tag22,,}_loci_conj_map.txt
        rm -f $snp2gene/${tag11,,}_${tag22,,}_loci_conj.txt
    fi
    if [ "${run_conj_flag:2:1}" = "Y" ]; then
        rm -f $snp2gene/${tag11,,}_${tag22,,}_genes_conj.txt
    fi
    if [ "${run_conj_flag:3:1}" = "Y" ]; then
        rm -f $gene2func/${tag11,,}_${tag22,,}_go_conj.txt
    fi
    if [ "${run_conj_flag:4:1}" = "Y" ]; then
        rm -f $gene2func/${tag11,,}_${tag22,,}_path_conj.txt
    fi
fi

if [ -z "${run_cond_flag##*Y*}" ]; then
    snp2gene2=$fuma_outfolder/snp2gene2
    mkdir -p $snp2gene2
    if [ ${#run_cond_flag} -eq 2 ]; then
        if [ "${run_cond_flag:0:1}" = "Y" ]; then
            rm -f $snp2gene2/${tag11,,}_${tag22,,}_snps_cond.txt
            rm -f $snp2gene2/${tag22,,}_${tag11,,}_snps_cond.txt
        fi
        if [ "${run_cond_flag:1:1}" = "Y" ]; then
            rm -f $snp2gene2/${tag11,,}_${tag22,,}_loci_cond_map_2.txt
            rm -f $snp2gene2/${tag22,,}_${tag11,,}_loci_cond_map_2.txt
            rm -f $snp2gene2/${tag11,,}_${tag22,,}_loci_cond_2.txt
            rm -f $snp2gene2/${tag22,,}_${tag11,,}_loci_cond_2.txt
        fi
    fi

    if [ ${#run_cond_flag} -eq 5 ]; then
        gene2func2=$fuma_outfolder/gene2func2
        mkdir -p $gene2func2
        if [ "${run_cond_flag:0:1}" = "Y" ]; then
            rm -f $snp2gene2/${tag11,,}_${tag22,,}_snps_cond.txt
            rm -f $snp2gene2/${tag22,,}_${tag11,,}_snps_cond.txt
        fi
        if [ "${run_cond_flag:1:1}" = "Y" ]; then
            rm -f $snp2gene2/${tag11,,}_${tag22,,}_loci_cond_map.txt
            rm -f $snp2gene2/${tag22,,}_${tag11,,}_loci_cond_map.txt
            rm -f $snp2gene2/${tag11,,}_${tag22,,}_loci_cond.txt
            rm -f $snp2gene2/${tag22,,}_${tag11,,}_loci_cond.txt
        fi
        if [ "${run_cond_flag:2:1}" = "Y" ]; then
            rm -f $snp2gene2/${tag11,,}_${tag22,,}_genes_cond.txt
            rm -f $snp2gene2/${tag22,,}_${tag11,,}_genes_cond.txt
        fi
        if [ "${run_cond_flag:3:1}" = "Y" ]; then
            rm -f $gene2func2/${tag11,,}_${tag22,,}_go_cond.txt
            rm -f $gene2func2/${tag22,,}_${tag11,,}_go_cond.txt
        fi
        if [ "${run_cond_flag:4:1}" = "Y" ]; then
            rm -f $gene2func2/${tag11,,}_${tag22,,}_path_cond.txt
            rm -f $gene2func2/${tag22,,}_${tag11,,}_path_cond.txt
        fi
    fi
fi

run_flag='Y'
n=1
for tag in `echo $taglist | awk -F ':' '{print $2}' | sed 's/,/\n/g'`; do
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

    n=$((n+1))

    echo "---------------------------------------------------------"
    echo $tag1 '#' $tag2
    echo $sm1 '#' $sm2
    echo $keyword1 '#' $keyword2
    echo $oldloci1 '#' $oldloci2
    echo "---------------------------------------------------------"

    if [ "$oldloci1" != "-" ] && [ ! -f $oldloci1 ]; then
        echo "NOTE: $oldloci1 doesn't exist"
        run_flag='N'
        continue
    fi

    if [ "$oldloci2" != "-" ] && [ ! -f $oldloci2 ]; then
        echo "NOTE: $oldloci2 doesn't exist"
        run_flag='N'
        continue
    fi

    #-------------------------------conjFDR-------------------------------#
    if [ -z "${run_conj_flag##*Y*}" ]; then
        conjfdr=$pleiofdr_outfolder/${sm1}_vs_${sm2}_conjfdr
        if [ ! -d $conjfdr ]; then
            echo "NOTE: $conjfdr doesn't exist, try to switch rev_flag"
            run_flag='N'
            continue
        fi

        fuma_conj=$snp2gene/${tag1}_${tag2}_conj
        fuma_conj2=$gene2func/${tag1}_${tag2}_conj

        if [ ! -d $fuma_conj ]; then
            echo "NOTE: $fuma_conj doesn't exist, please place fuma conjFDR snp2gene results into it"
            run_flag='N'
            continue
        fi

        if [ ! -d $fuma_conj2 ]; then
            echo "NOTE: $fuma_conj2 doesn't exist, please place fuma conjFDR gene2func results into it"
            run_flag='N'
            continue
        fi
    fi

    if [ -z "${run_cond_flag##*Y*}" ]; then
        condfdr=$pleiofdr_outfolder/${sm1}_vs_${sm2}_condfdr
        if [ ! -d $condfdr ]; then
            echo "NOTE: $condfdr doesn't exist"
            run_flag='N'
            continue
        fi

        condfdr2=$pleiofdr_outfolder/${sm2}_vs_${sm1}_condfdr
        if [ ! -d $condfdr2 ]; then
            echo "NOTE: $condfdr2 doesn't exist"
            run_flag='N'
            continue
        fi

        if [ ${#run_cond_flag} -eq 2 ]; then
            fuma_cond=$snp2gene2/${tag11}_${tag22}_cond
            if [ ${#run_cond_flag} -eq 2 ] && [ ! -d $fuma_cond ]; then
                echo "NOTE: $fuma_cond doesn't exist, please place fuma condFDR snp2gene results into it"
                exit 1
            fi
        fi

        if [ ${#run_cond_flag} -eq 5 ]; then
            fuma_cond_11=$snp2gene2/${tag1}_${tag2}_cond
            fuma_cond_12=$snp2gene2/${tag2}_${tag1}_cond
            fuma_cond_21=$gene2func2/${tag1}_${tag2}_cond
            fuma_cond_22=$gene2func2/${tag2}_${tag1}_cond

            if [ ! -d $fuma_cond_11 ]; then
                echo "NOTE: $fuma_cond_11 doesn't exist, please place fuma condFDR snp2gene results into it"
                run_flag='N'
                continue
            fi

            if [ ! -d $fuma_cond_12 ]; then
                echo "NOTE: $fuma_cond_12 doesn't exist, please place fuma condFDR snp2gene results into it"
                run_flag='N'
                continue
            fi

            if [ ! -d $fuma_cond_21 ]; then
                echo "NOTE: $fuma_cond_21 doesn't exist, please place fuma condFDR gene2func results into it"
                run_flag='N'
                continue
            fi

            if [ ! -d $fuma_cond_22 ]; then
                echo "NOTE: $fuma_cond_22 doesn't exist, please place fuma condFDR gene2func results into it"
                run_flag='N'
                continue
            fi
        fi
    fi

    if [ -z "${run_conj_flag##*Y*}" ]; then
        if [ "${run_conj_flag:0:1}" = "Y" ]; then
            echo "compiling conj snp table snp2gene/${tag1}_vs_${tag2}_snps.txt ..."
            sh $(dirname $0)/../tables/fdr_fuma_snp_table.sh $conjfdr/conj.result.clump.snps.csv 0.1 0.6 $fuma_conj/snps.txt $sumstat_folder/${sm1}.sumstats.gz $sumstat_folder/${sm2}.sumstats.gz $tag1 $tag2 $snp2gene
            echo "$tag1 & $tag2" >> $snp2gene/${tag11,,}_${tag22,,}_snps_conj.txt
            cat $snp2gene/${tag1}_vs_${tag2}_snps.txt >> $snp2gene/${tag11,,}_${tag22,,}_snps_conj.txt
            echo "" >> $snp2gene/${tag11,,}_${tag22,,}_snps_conj.txt
        fi
        if [ "${run_conj_flag:1:1}" = "Y" ]; then
            #gwasc=$fuma_conj/gwascatalog.txt
            echo "compiling conj loci table snp2gene/${tag1}_vs_${tag2}_loci.txt ..."
            sh $(dirname $0)/../tables/fdr_fuma_loci_table.sh $conjfdr/conj.result.clump.loci.csv $snp2gene/${tag1}_vs_${tag2}_snps.txt $snp2gene/${tag1}_vs_${tag2}_loci.txt $oldloci1 $oldloci2 $gwasc "$keyword1" "$keyword2"
            echo "$snp2gene/${tag1}_vs_${tag2}_loci.txt	$tag" >> $snp2gene/${tag11,,}_${tag22,,}_loci_conj_map.txt
            echo "$tag1 & $tag2" >> $snp2gene/${tag11,,}_${tag22,,}_loci_conj.txt
            cat $snp2gene/${tag1}_vs_${tag2}_loci.txt >> $snp2gene/${tag11,,}_${tag22,,}_loci_conj.txt
            echo "" >> $snp2gene/${tag11,,}_${tag22,,}_loci_conj.txt
        fi
        if [ "${run_conj_flag:2:1}" = "Y" ]; then
            echo "compiling conj gene table snp2gene/${tag1}_vs_${tag2}_genes.txt ..."
            sh $(dirname $0)/../tables/fdr_fuma_gene_table.sh $snp2gene/${tag1}_vs_${tag2}_snps.txt $fuma_conj/genes.txt $snp2gene/${tag1}_vs_${tag2}_genes.txt
            echo "$tag1 & $tag2" >> snp2gene/${tag11,,}_${tag22,,}_genes_conj.txt
            cat $snp2gene/${tag1}_vs_${tag2}_genes.txt >> $snp2gene/${tag11,,}_${tag22,,}_genes_conj.txt
            echo "" >> $snp2gene/${tag11,,}_${tag22,,}_genes_conj.txt
        fi
        if [ "${run_conj_flag:3:1}" = "Y" ]; then
            echo "compiling conj gene ontology table gene2func/${tag1}_vs_${tag2}_go.txt ..."
            sh $(dirname $0)/../tables/fdr_fuma_go_table.sh $snp2gene/${tag1}_vs_${tag2}_genes.txt $fuma_conj2/GS.txt $gene2func/${tag1}_vs_${tag2}_go.txt
            echo "$tag1 & $tag2" >> $gene2func/${tag11,,}_${tag22,,}_go_conj.txt
            cat $gene2func/${tag1}_vs_${tag2}_go.txt >> $gene2func/${tag11,,}_${tag22,,}_go_conj.txt
            echo "" >> $gene2func/${tag11,,}_${tag22,,}_go_conj.txt
        fi
        if [ "${run_conj_flag:4:1}" = "Y" ]; then
            echo "compiling conj pathway table gene2func/${tag1}_vs_${tag2}_path.txt ..."
            sh $(dirname $0)/../tables/fdr_cpdb_pathway_table.sh $snp2gene/${tag1}_vs_${tag2}_genes.txt $gene2func/ORA_${tag1}_${tag2}.tab $gene2func/${tag1}_vs_${tag2}_path.txt
            echo "$tag1 & $tag2" >> $gene2func/${tag11,,}_${tag22,,}_path_conj.txt
            cat $gene2func/${tag1}_vs_${tag2}_path.txt >> $gene2func/${tag11,,}_${tag22,,}_path_conj.txt
            echo "" >> $gene2func/${tag11,,}_${tag22,,}_path_conj.txt
        fi
    fi
    #-------------------------------condFDR-------------------------------#
    if [ -z "${run_cond_flag##*Y*}" ]; then
        if [ ${#run_cond_flag} -eq 2 ]; then
            if [ "${run_cond_flag:0:1}" = "Y" ]; then
                echo "compiling cond snp table snp2gene2/${tag1}_vs_${tag2}_snps.txt ..."
                sh $(dirname $0)/../tables/fdr_fuma_snp_table.sh $condfdr/cond.result.clump.snps.csv 0.1 0.6 $fuma_cond/snps.txt $sumstat_folder/${sm1}.sumstats.gz $sumstat_folder/${sm2}.sumstats.gz $tag1 $tag2 $snp2gene2
                echo "compiling cond snp table snp2gene2/${tag2}_vs_${tag1}_snps.txt ..."
                sh $(dirname $0)/../tables/fdr_fuma_snp_table.sh $condfdr2/cond.result.clump.snps.csv 0.1 0.6 $fuma_cond/snps.txt $sumstat_folder/${sm2}.sumstats.gz $sumstat_folder/${sm1}.sumstats.gz $tag2 $tag1 $snp2gene2
            fi
            if [ "${run_cond_flag:1:1}" = "Y" ]; then
                #gwasc=$fuma_cond/gwascatalog.txt
                echo "compiling cond loci table snp2gene2/${tag1}_vs_${tag2}_loci.txt ..."
                sh $(dirname $0)/../tables/fdr_fuma_loci_table.sh $condfdr/cond.result.clump.loci.csv $snp2gene2/${tag1}_vs_${tag2}_snps.txt $snp2gene2/${tag1}_vs_${tag2}_loci_2.txt $oldloci1 $oldloci2 $gwasc "$keyword1" "$keyword2"
                echo "$snp2gene2/${tag1}_vs_${tag2}_loci_2.txt	$tag" >> $snp2gene2/${tag11,,}_${tag22,,}_loci_cond_map_2.txt
                echo "$tag1 | $tag2" >> $snp2gene2/${tag11,,}_${tag22,,}_loci_cond_2.txt
                cat $snp2gene2/${tag1}_vs_${tag2}_loci_2.txt >> $snp2gene2/${tag11,,}_${tag22,,}_loci_cond_2.txt
                echo "" >> $snp2gene2/${tag11,,}_${tag22,,}_loci_cond_2.txt

                echo "compiling cond loci table snp2gene2/${tag2}_vs_${tag1}_loci.txt ..."
                sh $(dirname $0)/../tables/fdr_fuma_loci_table.sh $condfdr/cond.result.clump.loci.csv $snp2gene2/${tag2}_vs_${tag1}_snps.txt $snp2gene2/${tag2}_vs_${tag1}_loci_2.txt $oldloci2 $oldloci1 $gwasc "$keyword2" "$keyword1"
                echo "$snp2gene2/${tag2}_vs_${tag1}_loci_2.txt	$tag" >> $snp2gene2/${tag22,,}_${tag11,,}_loci_cond_map_2.txt
                echo "$tag2 | $tag1" >> $snp2gene2/${tag22,,}_${tag11,,}_loci_cond_2.txt
                cat $snp2gene2/${tag2}_vs_${tag1}_loci_2.txt >> $snp2gene2/${tag22,,}_${tag11,,}_loci_cond_2.txt
                echo "" >> $snp2gene2/${tag22,,}_${tag11,,}_loci_cond_2.txt
            fi
        fi

        if [ ${#run_cond_flag} -eq 5 ]; then
            if [ "${run_cond_flag:0:1}" = "Y" ]; then
                echo "compiling cond snp table snp2gene2/${tag1}_vs_${tag2}_snps.txt ..."
                sh $(dirname $0)/../tables/fdr_fuma_snp_table.sh $condfdr/cond.result.clump.snps.csv 0.1 0.6 $fuma_cond_11/snps.txt $sumstat_folder/${sm1}.sumstats.gz $sumstat_folder/${sm2}.sumstats.gz $tag1 $tag2 $snp2gene2
                echo "$tag1 | $tag2" >> $snp2gene2/${tag11,,}_${tag22,,}_snps_cond.txt
                cat $snp2gene2/${tag1}_vs_${tag2}_snps.txt >> $snp2gene2/${tag11,,}_${tag22,,}_snps_cond.txt
                echo "" >> $snp2gene2/${tag11,,}_${tag22,,}_snps_cond.txt

                echo "compiling cond snp table snp2gene2/${tag2}_vs_${tag1}_snps.txt ..."
                sh $(dirname $0)/../tables/fdr_fuma_snp_table.sh $condfdr/cond.result.clump.snps.csv 0.1 0.6 $fuma_cond_12/snps.txt $sumstat_folder/${sm2}.sumstats.gz $sumstat_folder/${sm1}.sumstats.gz $tag2 $tag1 $snp2gene2
                echo "$tag2 | $tag1" >> $snp2gene2/${tag22,,}_${tag11,,}_snps_cond.txt
                cat $snp2gene2/${tag2}_vs_${tag1}_snps.txt >> $snp2gene2/${tag22,,}_${tag11,,}_snps_cond.txt
                echo "" >> $snp2gene2/${tag22,,}_${tag11,,}_snps_cond.txt
            fi
            if [ "${run_cond_flag:1:1}" = "Y" ]; then
                #gwasc=$fuma_cond_11/gwascatalog.txt
                echo "compiling cond loci table snp2gene2/${tag1}_vs_${tag2}_loci.txt ..."
                sh $(dirname $0)/../tables/fdr_fuma_loci_table.sh $condfdr/cond.result.clump.loci.csv $snp2gene2/${tag1}_vs_${tag2}_snps.txt $snp2gene2/${tag1}_vs_${tag2}_loci.txt $oldloci1 $oldloci2 $gwasc "$keyword1" "$keyword2"
                echo "$snp2gene2/${tag1}_vs_${tag2}_loci.txt	$tag" >> $snp2gene2/${tag11,,}_${tag22,,}_loci_cond_map.txt
                echo "$tag1 | $tag2" >> $snp2gene2/${tag11,,}_${tag22,,}_loci_cond.txt
                cat $snp2gene2/${tag1}_vs_${tag2}_loci.txt >> $snp2gene2/${tag11,,}_${tag22,,}_loci_cond.txt
                echo "" >> $snp2gene2/${tag11,,}_${tag22,,}_loci_cond.txt

                #gwasc=$fuma_cond_12/gwascatalog.txt
                echo "compiling cond loci table snp2gene2/${tag2}_vs_${tag1}_loci.txt ..."
                sh $(dirname $0)/../tables/fdr_fuma_loci_table.sh $condfdr/cond.result.clump.loci.csv $snp2gene2/${tag2}_vs_${tag1}_snps.txt $snp2gene2/${tag2}_vs_${tag1}_loci.txt $oldloci2 $oldloci1 $gwasc "$keyword2" "$keyword1"
                echo "$snp2gene2/${tag2}_vs_${tag1}_loci.txt	$tag" >> $snp2gene2/${tag22,,}_${tag11,,}_loci_cond_map.txt
                echo "$tag2 | $tag1" >> $snp2gene2/${tag22,,}_${tag11,,}_loci_cond.txt
                cat $snp2gene2/${tag2}_vs_${tag1}_loci.txt >> $snp2gene2/${tag22,,}_${tag11,,}_loci_cond.txt
                echo "" >> $snp2gene2/${tag22,,}_${tag11,,}_loci_cond.txt
            fi
            if [ "${run_cond_flag:2:1}" = "Y" ]; then
                echo "compiling cond gene table snp2gene2/${tag1}_vs_${tag2}_genes.txt ..."
                sh $(dirname $0)/../tables/fdr_fuma_gene_table.sh $snp2gene2/${tag1}_vs_${tag2}_snps.txt $fuma_cond_11/genes.txt $snp2gene2/${tag1}_vs_${tag2}_genes.txt
                echo "$tag1 | $tag2" >> $snp2gene2/${tag11,,}_${tag22,,}_genes_cond.txt
                cat $snp2gene2/${tag1}_vs_${tag2}_genes.txt >> $snp2gene2/${tag11,,}_${tag22,,}_genes_cond.txt
                echo "" >> $snp2gene2/${tag11,,}_${tag22,,}_genes_cond.txt

                echo "compiling cond gene table snp2gene2/${tag2}_vs_${tag1}_genes.txt ..."
                sh $(dirname $0)/../tables/fdr_fuma_gene_table.sh $snp2gene2/${tag2}_vs_${tag1}_snps.txt $fuma_cond_12/genes.txt $snp2gene2/${tag2}_vs_${tag1}_genes.txt
                echo "$tag2 | $tag1" >> $snp2gene2/${tag22,,}_${tag11,,}_genes_cond.txt
                cat $snp2gene2/${tag2}_vs_${tag1}_genes.txt >> $snp2gene2/${tag22,,}_${tag11,,}_genes_cond.txt
                echo "" >> $snp2gene2/${tag22,,}_${tag11,,}_genes_cond.txt
            fi
            if [ "${run_cond_flag:3:1}" = "Y" ]; then
                echo "compiling cond gene ontology table gene2func2/${tag1}_vs_${tag2}_go.txt ..."
                sh $(dirname $0)/../tables/fdr_fuma_go_table.sh $snp2gene2/${tag1}_vs_${tag2}_genes.txt $fuma_cond_21/GS.txt $gene2func2/${tag1}_vs_${tag2}_go.txt
                echo "$tag1 | $tag2" >> $gene2func2/${tag11,,}_${tag22,,}_go_cond.txt
                cat $gene2func2/${tag1}_vs_${tag2}_go.txt >> $gene2func2/${tag11,,}_${tag22,,}_go_cond.txt
                echo "" >> $gene2func2/${tag11,,}_${tag22,,}_go_cond.txt

                echo "compiling cond gene ontology table gene2func2/${tag2}_vs_${tag1}_go.txt ..."
                sh $(dirname $0)/../tables/fdr_fuma_go_table.sh $snp2gene2/${tag2}_vs_${tag1}_genes.txt $fuma_cond_22/GS.txt $gene2func2/${tag2}_vs_${tag1}_go.txt
                echo "$tag2 | $tag1" >> $gene2func2/${tag22,,}_${tag11,,}_go_cond.txt
                cat $gene2func2/${tag2}_vs_${tag1}_go.txt >> $gene2func2/${tag22,,}_${tag11,,}_go_cond.txt
                echo "" >> $gene2func2/${tag22,,}_${tag11,,}_go_cond.txt
            fi
            if [ "${run_cond_flag:4:1}" = "Y" ]; then
                echo "compiling cond pathway table gene2func2/${tag1}_vs_${tag2}_path.txt ..."
                sh $(dirname $0)/../tables/fdr_cpdb_pathway_table.sh $snp2gene2/${tag1}_vs_${tag2}_genes.txt $gene2func2/ORA_${tag1}_${tag2}.tab $gene2func2/${tag1}_vs_${tag2}_path.txt
                echo "$tag1 | $tag2" >> $gene2func2/${tag11,,}_${tag22,,}_path_cond.txt
                cat $gene2func2/${tag1}_vs_${tag2}_path.txt >> $gene2func2/${tag11,,}_${tag22,,}_path_cond.txt
                echo "" >> $gene2func2/${tag11,,}_${tag22,,}_path_cond.txt

                echo "compiling cond pathway table gene2func2/${tag2}_vs_${tag1}_path.txt ..."
                sh $(dirname $0)/../tables/fdr_cpdb_pathway_table.sh $snp2gene2/${tag2}_vs_${tag1}_genes.txt $gene2func2/ORA_${tag2}_${tag1}.tab $gene2func2/${tag2}_vs_${tag1}_path.txt
                echo "$tag2 | $tag1" >> $gene2func2/${tag22,,}_${tag11,,}_path_cond.txt
                cat $gene2func2/${tag2}_vs_${tag1}_path.txt >> $gene2func2/${tag22,,}_${tag11,,}_path_cond.txt
                echo "" >> $gene2func2/${tag22,,}_${tag11,,}_path_cond.txt
            fi
        fi
    fi
    #---------------------------------------------------------------------#
done

if [ "$run_flag" = "N" ]; then
    exit 0
fi

#-------------------------------conjFDR-------------------------------#
if [ "${run_conj_flag:1:1}" = "Y" ]; then
    if [ -s $snp2gene/${tag11,,}_${tag22,,}_loci_conj_map.txt ] && [ `cat $snp2gene/${tag11,,}_${tag22,,}_loci_conj_map.txt | wc -l` -eq $((n_tag-2)) ]; then
        sh $(dirname $0)/../tables/fdr_fuma_loci_overlap.sh $snp2gene/${tag11,,}_${tag22,,}_loci_conj_map.txt $snp2gene/${tag11,,}_${tag22,,}_conj_shared_loci.txt
        rm -f $snp2gene/${tag11,,}_${tag22,,}_loci_conj.txt
        for tag in `echo $taglist | awk -F ':' '{print $2}' | sed 's/,/\n/g'`; do
            if [ "$rev_flag" = "Y" ]; then
                tag1=$tag
                tag2=$tagA
            else
                tag1=$tagA
                tag2=$tag
            fi
            echo "$tag1 & $tag2" >> $snp2gene/${tag11,,}_${tag22,,}_loci_conj.txt
            cat $snp2gene/${tag1}_vs_${tag2}_loci.txt >> $snp2gene/${tag11,,}_${tag22,,}_loci_conj.txt
            echo "" >> $snp2gene/${tag11,,}_${tag22,,}_loci_conj.txt
        done
    fi
fi
#-------------------------------condFDR-------------------------------#
if [ "${run_cond_flag:1:1}" = "Y" ]; then
    if [ -s $snp2gene2/${tag11,,}_${tag22,,}_loci_cond_map.txt ] && [ `cat $snp2gene2/${tag11,,}_${tag22,,}_loci_cond_map.txt | wc -l` -eq $((n_tag-2)) ]; then
        sh $(dirname $0)/../tables/fdr_fuma_loci_overlap.sh $snp2gene2/${tag11,,}_${tag22,,}_loci_cond_map.txt $snp2gene2/${tag11,,}_${tag22,,}_cond_shared_loci.txt
        rm -f $snp2gene2/${tag11,,}_${tag22,,}_loci_cond.txt
        for tag in `echo $taglist | awk -F ':' '{print $2}' | sed 's/,/\n/g'`; do
            if [ "$rev_flag" = "Y" ]; then
                tag1=$tag
                tag2=$tagA
            else
                tag1=$tagA
                tag2=$tag
            fi
            echo "$tag1 | $tag2" >> $snp2gene2/${tag11,,}_${tag22,,}_loci_cond.txt
            cat $snp2gene2/${tag1}_vs_${tag2}_loci.txt >> $snp2gene2/${tag11,,}_${tag22,,}_loci_cond.txt
            echo "" >> $snp2gene2/${tag11,,}_${tag22,,}_loci_cond.txt
        done
    fi
    if [ -s $snp2gene2/${tag22,,}_${tag11,,}_loci_cond_map.txt ] && [ `cat $snp2gene2/${tag22,,}_${tag11,,}_loci_cond_map.txt | wc -l` -eq $((n_tag-2)) ]; then
        sh $(dirname $0)/../tables/fdr_fuma_loci_overlap.sh $snp2gene2/${tag22,,}_${tag11,,}_loci_cond_map.txt $snp2gene2/${tag22,,}_${tag11,,}_cond_shared_loci.txt
        rm -f $snp2gene2/${tag22,,}_${tag11,,}_loci_cond.txt
        for tag in `echo $taglist | awk -F ':' '{print $2}' | sed 's/,/\n/g'`; do
            if [ "$rev_flag" = "Y" ]; then
                tag1=$tag
                tag2=$tagA
            else
                tag1=$tagA
                tag2=$tag
            fi
            echo "$tag2 | $tag1" >> $snp2gene2/${tag22,,}_${tag11,,}_loci_cond.txt
            cat $snp2gene2/${tag2}_vs_${tag1}_loci.txt >> $snp2gene2/${tag22,,}_${tag11,,}_loci_cond.txt
            echo "" >> $snp2gene2/${tag22,,}_${tag11,,}_loci_cond.txt
        done
    fi
    if [ -s $snp2gene2/${tag11,,}_${tag22,,}_loci_cond_map_2.txt ] && [ `cat $snp2gene2/${tag11,,}_${tag22,,}_loci_cond_map_2.txt | wc -l` -eq $((n_tag-2)) ]; then
        sh $(dirname $0)/../tables/fdr_fuma_loci_overlap.sh $snp2gene2/${tag11,,}_${tag22,,}_loci_cond_map_2.txt $snp2gene2/${tag11,,}_${tag22,,}_cond_shared_loci_2.txt
        rm -f $snp2gene2/${tag11,,}_${tag22,,}_loci_cond_2.txt
        for tag in `echo $taglist | awk -F ':' '{print $2}' | sed 's/,/\n/g'`; do
            if [ "$rev_flag" = "Y" ]; then
                tag1=$tag
                tag2=$tagA
            else
                tag1=$tagA
                tag2=$tag
            fi
            echo "$tag1 | $tag2" >> $snp2gene2/${tag11,,}_${tag22,,}_loci_cond_2.txt
            cat $snp2gene2/${tag1}_vs_${tag2}_loci_2.txt >> $snp2gene2/${tag11,,}_${tag22,,}_loci_cond_2.txt
            echo "" >> $snp2gene2/${tag11,,}_${tag22,,}_loci_cond_2.txt
        done
    fi
    if [ -s $snp2gene2/${tag22,,}_${tag11,,}_loci_cond_map_2.txt ] && [ `cat $snp2gene2/${tag22,,}_${tag11,,}_loci_cond_map_2.txt | wc -l` -eq $((n_tag-2)) ]; then
        sh $(dirname $0)/../tables/fdr_fuma_loci_overlap.sh $snp2gene2/${tag22,,}_${tag11,,}_loci_cond_map_2.txt $snp2gene2/${tag22,,}_${tag11,,}_cond_shared_loci_2.txt
        rm -f $snp2gene2/${tag22,,}_${tag11,,}_loci_cond_2.txt
        for tag in `echo $taglist | awk -F ':' '{print $2}' | sed 's/,/\n/g'`; do
            if [ "$rev_flag" = "Y" ]; then
                tag1=$tag
                tag2=$tagA
            else
                tag1=$tagA
                tag2=$tag
            fi
            echo "$tag2 | $tag1" >> $snp2gene2/${tag22,,}_${tag11,,}_loci_cond_2.txt
            cat $snp2gene2/${tag2}_vs_${tag1}_loci_2.txt >> $snp2gene2/${tag22,,}_${tag11,,}_loci_cond_2.txt
            echo "" >> $snp2gene2/${tag22,,}_${tag11,,}_loci_cond_2.txt
        done
    fi
fi
#---------------------------------------------------------------------#
