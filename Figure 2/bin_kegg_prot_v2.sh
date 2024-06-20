#!/usr/bin/sh
#A pipeline for bin kegg annotation
#jinh 20191119
echo "V2.0, last updata: 20191119, jinhao94@126.com"
ext=faa
echo -e "Default extension is \033[31m$ext \033[0m"
bin_dir=$1
out_dir=$2

if [ $# -eq 3 ];then
    ext=$3
fi

function helpmessage {
echo -e "Usage: bin_kegg_prot_v2.sh \033[33m bin_folder bin_KEGG_result extension\033[0m (\033[31mdefault is faa, the extension is very important!\033[0m)"
# echo -e "$0"
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


kegg_db=/nvmessdnode3/opt/database/KEGG_Archaea_Bacteria_fungi/KEGG_ABF.fasta.dmnd
mkdir -p $out_dir/protein_faa
mkdir $out_dir/combine_result
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~prodigal
# echo "Running prodigal"
# ls -d $bin_dir/*$ext | parallel --plus -j 32 prodigal -i {} -a $out_dir/protein_faa/{/.}.faa -o $out_dir/protein_faa/{/.}.gff -d $out_dir/protein_faa/{/.}.fnn -q -f gff 
# rename gene names and get gene length 
ls -d $bin_dir/*$ext | parallel -j 20 --plus 181120rename_prodigal_seq_names.py {} {/.} $out_dir/protein_faa/{/.}.rename_faa

#ls -d $bin_dir/*$ext | parallel -j 10 --plus cp {} $out_dir/protein_faa/{/.}.rename_faa
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~run diamond
echo "Running diamond..."
#cat files
cat $out_dir/protein_faa/*rename_faa > $out_dir/protein_faa/combine.faa
diamond blastp --threads 32  --db $kegg_db --query $out_dir/protein_faa/combine.faa --outfmt 6 qseqid sseqid stitle pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore --out $out_dir/combine_result/combine.faa.kegg_dia.result --masking 0 --query-cover 70 --id 30 --quiet

#we have a best hits 
less $out_dir/combine_result/combine.faa.kegg_dia.result | perl -ne '@s=split /\s+/;print $_ unless $a eq $s[0];$a=$s[0];' > $out_dir/combine_result/combine.faa.kegg_dia.result.best_hits

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~get ko
echo "Getting ko..............."
perl -e 'open DB, "/nvmessdnode3/opt/database/KEGG_Archaea_Bacteria_fungi/Archaea_Bacteria_Fungi.kegg.ko.list"; while(<DB>){chomp; @s=split /\s+/; $h{@s[0]}=@s[1]}; open IN, @ARGV[0]; while(<IN>){chomp; @s=split /\t/; @s[1]=~/.*:(.*)/; print "$_\t$h{@s[1]}\n" if exists $h{@s[1]}}' $out_dir/combine_result/combine.faa.kegg_dia.result.best_hits > $out_dir/combine_result/combine.faa.kegg_dia.result.best_hits.ko

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~split result
echo "Splitting result.............."
mkdir $out_dir/Split_result
mkdir ${out_dir}/gene_info
perl -e 'open IN1, @ARGV[0]; while(<IN1>){chomp; @s=split /\s+/; @s1 = split /_/, @s[0]; $a=(join "_", @s1[0..$#s1]); $o=@ARGV[1]."/".$a; if($a ne $n){open O, ">$o"; print O "$_\n" ; $n=$a}else{print O "$_\n"}}' $out_dir/combine_result/combine.faa.kegg_dia.result.best_hits.ko $out_dir/Split_result

ls -d $out_dir/Split_result/* | parallel  --plus -j 3 "wc -l {} > $out_dir/gene_info/{/}.sm; cut -f16 {} | sort -u | wc -l >> $out_dir/gene_info/{/}.sm ; grep '^>' $out_dir/protein_faa/{/}.rename_faa -c >> $out_dir/gene_info/{/}.sm"
# mv $out_dir/Split_result/*.sm $out_dir/gene_info

ls -d $out_dir/Split_result/* | parallel  -j 10 "cut -f16 {} | sort | uniq -c | awk '{print \$2\"\t\"\$1}' | sort -k2 | sed '1i KO\t{/}' > {}.gene "

# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~get ko matrix
echo "Getting ko matrix..."
cat $out_dir/Split_result/*gene | cut -f1 |sort -u > $out_dir/KO_raw_names

perl -e '%h; open IN, @ARGV[0]; while(<IN>){ chomp $_; $h{$_}=(); } ; opendir DI, @ARGV[1] ; foreach $file(readdir DI){next unless $file=~/gene/; $infile=@ARGV[1]."/".$file; open IN, $infile ; undef %h2; while(<IN>){chomp; @s=split /\t/; $h2{$s[0]}=$s[1]}; foreach $k(keys %h){if(exists $h2{$k}){push@{$h{$k}},$h2{$k};}else{ push@{$h{$k}},0 } } };foreach $k2(keys %h){@l=@{$h{$k2}};$o=(join "\t",@l); print "$k2\t$o\n"}' $out_dir/KO_raw_names $out_dir/Split_result | sort -k1 -r > $out_dir/bin_kegg_matrix

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~get bin annotation summary
for i in `ls -d ${out_dir}/gene_info/*sm` ; do perl -e 'open IN, @ARGV[0]; while(<IN>){chomp;@s=split /\s+/; push @p,$s[0]}; $o=join("\t",@p); print "@ARGV[0]\t$o\n"' $i > ${i}.f ; done
echo -e "bin\tnum_seq_get_ko\tnum_ko_category\tnum_protein" > $out_dir/bin_gene_info_summary 
ls -d $out_dir/gene_info/*sm.f | parallel cat {} '>>' $out_dir/bin_gene_info_summary
sed  -i 's/.\/\/gene_info\/\|\.sm//g' $out_dir/bin_gene_info_summary
echo "pipeline done..............."
