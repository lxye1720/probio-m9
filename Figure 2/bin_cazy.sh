#!/usr/bin/sh
#A pipeline for bin kegg annotation
#jinh 2019.5.12
ext=fa
dbdir=/nvmessdnode3/opt/database/CAZY_db
echo "Default extension is $ext"
bin_dir=$1
out_dir=$2

if [ $# -eq 3 ];then
    ext=$3
fi

function helpmessage {
echo "Usage bin_cazy.sh bin_folder bin_CAZY_result fa"
}
if [ $# -eq 0 ];then
    echo "Please input parameters!"
    helpmessage
    exit
fi

if [ ! -d $bin_dir ];then
    echo "bin folder not exist, exit..."
    helpmessage
    exit
fi

if [  -d $out_dir ];then
    echo "output folder exist, so remove..."
    rm -rf $out_dir
fi


echo "Running hmmer"
mkdir -p ${out_dir}/tmp_dir ${out_dir}/combine_result
ls -d ${bin_dir}/*.${ext} | parallel -j 8 python /nvmessdnode3/opt/links/run_dbcan.py2 --hmm_cpu 6 --db_dir /nvmessdnode3/opt/database/CAZY_db -t hmmer --out_dir ${out_dir}/tmp_dir/{/.} {} prok 

ls -d ${out_dir}/tmp_dir/* | parallel "sed -i 's/.hmm//g' {}/hmmer.out"
#best hits
ls  -d ${out_dir}/tmp_dir/* | parallel -j 8 "sed -i '1d' {}/hmmer.out"
ls  -d ${out_dir}/tmp_dir/* | parallel -j 8 "sort -k3 {}/hmmer.out > {}/hmmer.out.s "
ls  -d ${out_dir}/tmp_dir/* | parallel -j 8 "less {}/hmmer.out.s | perl -ne 'chomp; @s=split /\t/; if(\$s[0] ne \$n){print \"\$_\\n\"; \$n=\$s[0]; next  }' > {}/hmmer.out "

ls -d ${out_dir}/tmp_dir/* | parallel --plus "nm={}; onm=\${nm##*/} ; cut -f 1 {}/hmmer.out | sort | uniq -c | awk -F\" \" '{print \$2\"\t\"\$1}' > {}/\${onm}.hmm_prof"

rawpath=$PWD
cd ${out_dir}
ls -d tmp_dir/* | parallel mv {}/*.hmm_prof combine_result


#get header
ls -d combine_result/* | parallel 'sed -i "1 i111111Enzyme_id\t{/.}" {} '
cat combine_result/*.hmm_prof > all_com
sort all_com |  cut -f1 | uniq > all_com.uniq.names

less all_com.uniq.names | perl -e '%h;while(<>){ chomp $_; $h{$_}=(); } ; opendir DI,"./combine_result" ; foreach $file(readdir DI){next unless $file=~/hmm_prof/; $infile="./combine_result/".$file; open IN, $infile ; undef %h2; ;while(<IN>){chomp; @s=split /\t/; $h2{@s[0]}=@s[1]}; foreach $k(keys %h){if(exists $h2{$k}){push@{$h{$k}},$h2{$k};}else{ push@{$h{$k}},0 } } };foreach $k2(keys %h){@l=@{$h{$k2}};$o=(join "\t",@l); print "$k2\t$o\n"}' | sort -k1 | sed 's/111111Enzyme_id/Enzyme_id/g' > bin_cazy_result

#get info
echo "get info"
mkdir count_files
ls -d combine_result/*.hmm_prof | parallel echo {/.} '>' count_files/{/.}.me_count
ls -d tmp_dir/* | parallel 'sed "1d" {}/hmmer.out | wc -l >> count_files/{/}.me_count '
ls -d combine_result/*.hmm_prof | parallel "sed '1d' {} | wc -l >> count_files/{/.}.me_count"

for i in `ls -d count_files/*me_count ` ; do less $i | perl -e 'while(<>){chomp;@s=split /\s+/; push @p,$s[0]}; $o=join("\t",@p); print "$o\n"' > ${i}1 ; done

echo -e "bin\tget_enzyme_num\tenzyme_num" > bin_gene_info_summary 
ls -d count_files/*me_count1 | parallel cat {} '>>' bin_gene_info_summary
rm all_com all_com.uniq.names
cd $rawpath
#
echo "Done..................................."
