$HOSTNAME = ""
params.outdir = 'results'  

def pathChecker(input, path, type){
	cmd = "mkdir -p check && mv ${input} check/. "
	if (!input || input.empty()){
		input = file(path).getName().toString()
		cmd = "mkdir -p check && cd check && ln -s ${path} ${input} && cd .."
		if (path.indexOf('s3:') > -1 || path.indexOf('S3:') >-1){
			recursive = (type == "folder") ? "--recursive" : ""
			cmd = "mkdir -p check && cd check && aws s3 cp ${recursive} ${path} ${workDir}/${input} && ln -s ${workDir}/${input} . && cd .."
		} else if (path.indexOf('gs:') > -1 || path.indexOf('GS:') >-1){
			if (type == "folder"){
				cmd = "mkdir -p check ${workDir}/${input} && cd check && gsutil rsync -r ${path} ${workDir}/${input} && cp -R ${workDir}/${input} . && cd .."
			} else {
				cmd = "mkdir -p check && cd check && gsutil cp ${path} ${workDir}/${input} && cp -R ${workDir}/${input} . && cd .."
			}
		} else if (path.indexOf('/') == -1){
			cmd = ""
		}
}
	return [cmd,input]
}
if (!params.mate_split){params.mate_split = ""} 
if (!params.trans2gene_path){params.trans2gene_path = ""} 
if (!params.reads){params.reads = ""} 
if (!params.cellBarcodeFile){params.cellBarcodeFile = ""} 
if (!params.mate){params.mate = ""} 
// Stage empty file to be used as an optional input where required
ch_empty_file_1 = file("$baseDir/.emptyfiles/NO_FILE_1", hidden:true)
ch_empty_file_2 = file("$baseDir/.emptyfiles/NO_FILE_2", hidden:true)
ch_empty_file_3 = file("$baseDir/.emptyfiles/NO_FILE_3", hidden:true)
ch_empty_file_4 = file("$baseDir/.emptyfiles/NO_FILE_4", hidden:true)
ch_empty_file_5 = file("$baseDir/.emptyfiles/NO_FILE_5", hidden:true)

Channel.value(params.mate_split).into{g_8_1_g116_30;g_8_0_g116_31;g_8_1_g114_82;g_8_0_g114_131;g_8_2_g114_134;g_8_1_g186_11;g_8_1_g186_16;g_8_1_g186_21;g_8_1_g186_24;g_8_0_g186_28;g_8_0_g186_31;g_8_1_g186_23;g_8_1_g186_19;g_8_1_g186_20;g_8_1_g186_18}
g_125_0_g_161 = file(params.trans2gene_path, type: 'any')
g_125_0_g_175 = file(params.trans2gene_path, type: 'any')
if (params.reads){
Channel
	.fromFilePairs( params.reads , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 )
	.ifEmpty { error "Cannot find any fastq_set matching: ${params.reads}" }
	.set{g_158_1_g_166}
 } else {  
	g_158_1_g_166 = Channel.empty()
 }

g_159_0_g_166 = file(params.cellBarcodeFile, type: 'any')
Channel.value(params.mate).set{g_167_2_g_166}

//* params.gtf =  ""  //* @input
//* params.genome =  ""  //* @input
//* params.commondb =  ""  //* @input
//* params.genome_source =  ""  //* @input
//* params.gtf_source =  ""  //* @input
//* params.commondb_source =  ""  //* @input @optional

def downFile(path, task){
	println workDir
    if (path.take(1).indexOf("/") == 0){
      target=path
      if (task.executor == "awsbatch" || task.executor == "google-batch") {
      	a=file(path)
    	fname = a.getName().toString()
    	target = "${workDir}/${fname}"
    	if (!file(target).exists()){
    		a.copyTo(workDir)
    	}
      }
    } else {
      a=file(path)
      fname = a.getName().toString()
      target = "${workDir}/${fname}"
      if (!file(target).exists()){
    		a.copyTo(workDir)
      } 
    }
    return target
}

def getLastName (str){
	if (str.indexOf("/") > -1){
		return  str.substring(str.lastIndexOf('/')+1,str.length())
	} 
	return ""
}

process Check_and_Build_Module_Check_Genome_GTF {


output:
 file "${newNameFasta}"  into g184_21_genome00_g184_52, g184_21_genome01_g184_54
 file "${newNameGtf}"  into g184_21_gtfFile10_g184_53, g184_21_gtfFile10_g184_54

container 'quay.io/viascientific/pipeline_base_image:1.0'

when:
params.run_Download_Genomic_Sources == "yes"

script:
genomeSource = !file("${params.genome}").exists() ? params.genome_source : params.genome
genomeName = getLastName(genomeSource)

gtfSource = !file("${params.gtf}").exists() ? params.gtf_source : params.gtf
gtfName = getLastName(gtfSource)


newNameGtf = gtfName
newNameFasta = genomeName
if (gtfName.contains('.gz')) { newNameGtf =  newNameGtf - '.gz'  } 
if (genomeName.contains('.gz')) { newNameFasta =  newNameFasta - '.gz'  } 

runGzip = ""
if (gtfName.contains('.gz') || genomeName.contains('.gz')) {
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} 

slashCountGenome = params.genome_source.count("/")
cutDirGenome = slashCountGenome - 3;

slashCountGtf = params.gtf_source.count("/")
cutDirGtf = slashCountGtf - 3;

"""
if [ ! -e "${params.genome_source}" ] ; then
    echo "${params.genome_source} not found"
	if [[ "${params.genome_source}" =~ "s3" ]]; then
		echo "Downloading s3 path from ${params.genome_source}"
		aws s3 cp ${params.genome_source} ${workDir}/${genomeName} && ln -s ${workDir}/${genomeName} ${genomeName}
	elif [[ "${params.genome_source}" =~ "gs" ]]; then
		echo "Downloading gs path from ${params.genome_source}"
		gsutil cp  ${params.genome_source} ${workDir}/. && ln -s ${workDir}/${genomeName} ${genomeName}
	else
		echo "Downloading genome with wget"
		wget --no-check-certificate --secure-protocol=TLSv1 -l inf -nc -nH --cut-dirs=$cutDirGenome -R 'index.html*' -r --no-parent  ${params.genome_source}
	fi

else 
	ln -s ${params.genome_source} ${genomeName}
fi

if [ ! -e "${params.gtf_source}" ] ; then
    echo "${params.gtf_source} not found"
	if [[ "${params.gtf_source}" =~ "s3" ]]; then
		echo "Downloading s3 path from ${params.gtf_source}"
		aws s3 cp  ${params.gtf_source} ${workDir}/${gtfName} && ln -s ${workDir}/${gtfName} ${gtfName}
	elif [[ "${params.gtf_source}" =~ "gs" ]]; then
		echo "Downloading gs path from ${params.gtf_source}"
		gsutil cp  ${params.gtf_source} ${workDir}/. && ln -s ${workDir}/${gtfName} ${gtfName}
	else
		echo "Downloading gtf with wget"
		wget --no-check-certificate --secure-protocol=TLSv1 -l inf -nc -nH --cut-dirs=$cutDirGtf -R 'index.html*' -r --no-parent  ${params.gtf_source}
	fi

else 
	ln -s ${params.gtf_source} ${gtfName}
fi

$runGzip

"""




}

//* params.gtf2bed_path =  ""  //* @input
//* params.bed =  ""  //* @input

process Check_and_Build_Module_Check_BED12 {

input:
 file gtf from g184_21_gtfFile10_g184_53

output:
 file "${gtfName}.bed"  into g184_53_bed03_g184_54

when:
params.run_Download_Genomic_Sources == "yes"

script:
gtfName  = gtf.baseName
beddir = ""
if (params.bed.indexOf('/') > -1){
	beddir  = params.bed.substring(0, params.bed.lastIndexOf('/')) 
}
"""

if [ ! -e "${params.bed}" ] ; then
    echo "${params.bed} not found"
    perl ${params.gtf2bed_path} $gtf > ${gtfName}.bed
else 
	cp -n ${params.bed} ${gtfName}.bed
fi
if [ "${beddir}" != "" ] ; then
	mkdir -p ${beddir}
	cp -n ${gtfName}.bed ${params.bed} 
fi
"""




}

//* params.gtf2bed_path =  ""  //* @input
//* params.genome_sizes =  ""  //* @input

process Check_and_Build_Module_Check_chrom_sizes_and_index {

input:
 file genome from g184_21_genome00_g184_52

output:
 file "${genomeName}.chrom.sizes"  into g184_52_genomeSizes02_g184_54

when:
params.run_Download_Genomic_Sources == "yes"

script:
genomeName  = genome.baseName
genome_sizes_dir = ""
if (params.genome_sizes.indexOf('/') > -1){
	genome_sizes_dir  = params.genome_sizes.substring(0, params.genome_sizes.lastIndexOf('/')) 
}

"""
if [ ! -e "${params.genome_sizes}" ] ; then
    echo "${params.genome_sizes} not found"
    cat ${genome} | awk '\$0 ~ ">" {print c; c=0;printf substr(\$1,2,100) "\\t"; } \$0 !~ ">" {c+=length(\$0);} END { print c; }' > ${genomeName}.chrom.sizes
    ##clean first empty line
    sed -i '1{/^\$/d}' ${genomeName}.chrom.sizes
    if [ "${genome_sizes_dir}" != "" ] ; then
    	mkdir -p ${genome_sizes_dir}
		cp -n ${genomeName}.chrom.sizes ${params.genome_sizes} 
	fi
else 
	cp ${params.genome_sizes} ${genomeName}.chrom.sizes
fi

"""




}

g184_21_gtfFile10_g184_54= g184_21_gtfFile10_g184_54.ifEmpty([""]) 
g184_21_genome01_g184_54= g184_21_genome01_g184_54.ifEmpty([""]) 
g184_52_genomeSizes02_g184_54= g184_52_genomeSizes02_g184_54.ifEmpty([""]) 
g184_53_bed03_g184_54= g184_53_bed03_g184_54.ifEmpty([""]) 


process Check_and_Build_Module_check_files {

input:
 file gtf from g184_21_gtfFile10_g184_54
 file genome from g184_21_genome01_g184_54
 file genomeSizes from g184_52_genomeSizes02_g184_54
 file bed from g184_53_bed03_g184_54

output:
 file "*/${gtf2}" optional true  into g184_54_gtfFile01_g116_21
 file "*/${genome2}" optional true  into g184_54_genome10_g116_21
 file "*/${genomeSizes2}" optional true  into g184_54_genomeSizes22_g114_131, g184_54_genomeSizes21_g114_142
 file "*/${bed2}" optional true  into g184_54_bed31_g114_134

container 'quay.io/viascientific/pipeline_base_image:1.0'

script:
(cmd1, gtf2) = pathChecker(gtf, params.gtf, "file")
(cmd2, genome2) = pathChecker(genome, params.genome, "file")
(cmd3, genomeSizes2) = pathChecker(genomeSizes, params.genome_sizes, "file")
(cmd4, bed2) = pathChecker(bed, params.bed, "file")
"""
$cmd1
$cmd2
$cmd3
$cmd4
"""
}

build_STAR_index = params.STAR_Module_Check_Build_STAR_Index.build_STAR_index

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 5
    $MEMORY = 150
}
//* platform
//* platform
//* autofill

process STAR_Module_Check_Build_STAR_Index {

input:
 file genome from g184_54_genome10_g116_21
 file gtf from g184_54_gtfFile01_g116_21

output:
 file "STARIndex"  into g116_21_starIndex00_g116_26

when:
build_STAR_index == true && ((params.run_STAR && (params.run_STAR == "yes")) || !params.run_STAR)

script:
star_build_parameters = params.STAR_Module_Check_Build_STAR_Index.star_build_parameters
newDirName = "STARIndex" 
"""
if [ ! -e "${params.star_index}/SA" ] ; then
    echo "STAR index not found"
    mkdir -p $newDirName 
    STAR --runMode genomeGenerate ${star_build_parameters} --genomeDir $newDirName --genomeFastaFiles ${genome} --sjdbGTFfile ${gtf}
else 
	ln -s ${params.star_index} STARIndex
fi

"""





}

g116_21_starIndex00_g116_26= g116_21_starIndex00_g116_26.ifEmpty([""]) 

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 10
}
//* platform
//* platform
//* autofill
if (!((params.run_STAR && (params.run_STAR == "yes")) || !params.run_STAR)){
g116_21_starIndex00_g116_26.set{g116_26_starIndex02_g116_31}
} else {


process STAR_Module_check_STAR_files {

input:
 file star from g116_21_starIndex00_g116_26

output:
 file "*/${star2}" optional true  into g116_26_starIndex02_g116_31

container 'quay.io/viascientific/pipeline_base_image:1.0'

when:
(params.run_STAR && (params.run_STAR == "yes")) || !params.run_STAR

script:
(cmd, star2) = pathChecker(star, params.star_index, "folder")
"""
$cmd
"""
}
}


//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 4
}
//*

process extractValidReads_XL {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /validfastq\/.*$/) "validReads/$filename"}
input:
 file cellBarcode from g_159_0_g_166
 set val(name), file(reads) from g_158_1_g_166
 val mate from g_167_2_g_166

output:
 set val(name), file("validfastq/*")  into g_166_valid_fastq01_g186_28, g_166_valid_fastq00_g186_18

script:
readArr = reads.toString().split(" ")
"""
mkdir validfastq
path=\$(which extractValidReads_UMIspecial2.py) && cp \$path .
python extractValidReads_UMIspecial2.py -i ${readArr[0]} -o validfastq -b ${cellBarcode}
gzip validfastq/*
"""
}

//* params.run_Adapter_Removal =   "no"   //* @dropdown @options:"yes","no" @show_settings:"Adapter_Removal"
//* @style @multicolumn:{seed_mismatches, palindrome_clip_threshold, simple_clip_threshold} @condition:{Tool_for_Adapter_Removal="trimmomatic", seed_mismatches, palindrome_clip_threshold, simple_clip_threshold}, {Tool_for_Adapter_Removal="fastx_clipper", discard_non_clipped}


//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 5
    $MEMORY = 8
}
//* platform
//* platform
//* autofill
if (!((params.run_Adapter_Removal && (params.run_Adapter_Removal == "yes")) || !params.run_Adapter_Removal)){
g_166_valid_fastq00_g186_18.into{g186_18_reads01_g186_31; g186_18_reads00_g186_23}
g186_18_log_file10_g186_11 = Channel.empty()
} else {


process Adapter_Trimmer_Quality_Module_Adapter_Removal {

input:
 set val(name), file(reads) from g_166_valid_fastq00_g186_18
 val mate from g_8_1_g186_18

output:
 set val(name), file("reads/*.fastq.gz")  into g186_18_reads01_g186_31, g186_18_reads00_g186_23
 file "*.{fastx,trimmomatic}.log"  into g186_18_log_file10_g186_11

container 'quay.io/viascientific/trimmomatic:1.0'

when:
(params.run_Adapter_Removal && (params.run_Adapter_Removal == "yes")) || !params.run_Adapter_Removal

shell:
phred = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.phred
Tool_for_Adapter_Removal = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.Tool_for_Adapter_Removal
Adapter_Sequence = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.Adapter_Sequence
//trimmomatic_inputs
min_length = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.min_length
seed_mismatches = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.seed_mismatches
palindrome_clip_threshold = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.palindrome_clip_threshold
simple_clip_threshold = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.simple_clip_threshold

//fastx_clipper_inputs
discard_non_clipped = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.discard_non_clipped
    
remove_previous_reads = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.remove_previous_reads
discard_non_clipped_text = ""
if (discard_non_clipped == "yes") {discard_non_clipped_text = "-c"}
nameAll = reads.toString()
nameArray = nameAll.split(' ')
file2 = ""
if (nameAll.contains('.gz')) {
    newName =  nameArray[0] 
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1] }
} 
'''
#!/usr/bin/env perl
 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage;
 use Cwd qw();
 
runCmd("mkdir reads adapter unpaired");

open(OUT, ">adapter/adapter.fa");
my @adaps=split(/\n/,"!{Adapter_Sequence}");
my $i=1;
foreach my $adap (@adaps)
{
 print OUT ">adapter$i\\n$adap\\n";
 $i++;
}
close(OUT);

my $quality="!{phred}";
print "fastq quality: $quality\\n";
print "tool: !{Tool_for_Adapter_Removal}\\n";

if ("!{mate}" eq "pair") {
    if ("!{Tool_for_Adapter_Removal}" eq "trimmomatic") {
        runCmd("trimmomatic PE -threads !{task.cpus} -phred${quality} !{file1} !{file2} reads/!{name}.1.fastq.gz unpaired/!{name}.1.fastq.unpaired.gz reads/!{name}.2.fastq.gz unpaired/!{name}.2.fastq.unpaired.gz ILLUMINACLIP:adapter/adapter.fa:!{seed_mismatches}:!{palindrome_clip_threshold}:!{simple_clip_threshold} MINLEN:!{min_length} 2> !{name}.trimmomatic.log");
    } elsif ("!{Tool_for_Adapter_Removal}" eq "fastx_clipper") {
        print "Fastx_clipper is not suitable for paired reads.";
    }
} else {
    if ("!{Tool_for_Adapter_Removal}" eq "trimmomatic") {
        runCmd("trimmomatic SE -threads !{task.cpus}  -phred${quality} !{file1} reads/!{name}.fastq.gz ILLUMINACLIP:adapter/adapter.fa:!{seed_mismatches}:!{palindrome_clip_threshold}:!{simple_clip_threshold} MINLEN:!{min_length} 2> !{name}.trimmomatic.log");
    } elsif ("!{Tool_for_Adapter_Removal}" eq "fastx_clipper") {
        runCmd("fastx_clipper  -Q $quality -a !{Adapter_Sequence} -l !{min_length} !{discard_non_clipped_text} -v -i !{file1} -o reads/!{name}.fastq.gz > !{name}.fastx.log");
    }
}
if ("!{remove_previous_reads}" eq "true") {
    my $currpath = Cwd::cwd();
    my @paths = (split '/', $currpath);
    splice(@paths, -2);
    my $workdir= join '/', @paths;
    splice(@paths, -1);
    my $inputsdir = join '/', @paths;
    $inputsdir .= "/work";
    print "INFO: inputs reads will be removed if they are located in the $workdir $inputsdir\\n";
    my @listOfFiles = `readlink -e !{file1} !{file2}`;
    foreach my $targetFile (@listOfFiles){
        if (index($targetFile, $workdir) != -1 || index($targetFile, $inputsdir) != -1) {
            runCmd("rm -f $targetFile");
            print "INFO: $targetFile deleted.\\n";
        }
    }
}


##Subroutines
sub runCmd {
    my ($com) = @_;
    if ($com eq ""){
		return "";
    }
    my $error = system(@_);
    if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}
'''

}
}


//* @style @condition:{single_or_paired_end_reads="single", barcode_pattern1,remove_duplicates_based_on_UMI}, {single_or_paired_end_reads="pair", barcode_pattern1,barcode_pattern2}

if (!(params.run_UMIextract == "yes")){
g186_18_reads00_g186_23.set{g186_23_reads00_g186_19}
g186_23_log_file10_g186_24 = Channel.empty()
} else {


process Adapter_Trimmer_Quality_Module_UMIextract {

input:
 set val(name), file(reads) from g186_18_reads00_g186_23
 val mate from g_8_1_g186_23

output:
 set val(name), file("result/*.fastq.gz")  into g186_23_reads00_g186_19
 file "${name}.*.log"  into g186_23_log_file10_g186_24

label 'fastq_preprocessing'

when:
params.run_UMIextract == "yes" 

script:
readArray = reads.toString().split(' ')
file2 = ""
file1 =  readArray[0]
if (mate == "pair") {file2 =  readArray[1]}


single_or_paired_end_reads = params.Adapter_Trimmer_Quality_Module_UMIextract.single_or_paired_end_reads
barcode_pattern1 = params.Adapter_Trimmer_Quality_Module_UMIextract.barcode_pattern1
barcode_pattern2 = params.Adapter_Trimmer_Quality_Module_UMIextract.barcode_pattern2
UMIqualityFilterThreshold = params.Adapter_Trimmer_Quality_Module_UMIextract.UMIqualityFilterThreshold
phred = params.Adapter_Trimmer_Quality_Module_UMIextract.phred
remove_duplicates_based_on_UMI = params.Adapter_Trimmer_Quality_Module_UMIextract.remove_duplicates_based_on_UMI

"""
set +e
source activate umi_tools_env 2> /dev/null || true
mkdir result
if [ "${mate}" == "pair" ]; then
umi_tools extract --bc-pattern='${barcode_pattern1}' \
                  --bc-pattern2='${barcode_pattern2}' \
                  --extract-method=regex \
                  --stdin=${file1} \
                  --stdout=result/${name}_R1.fastq.gz \
                  --read2-in=${file2} \
                  --read2-out=result/${name}_R2.fastq.gz\
				  --quality-filter-threshold=${UMIqualityFilterThreshold} \
				  --quality-encoding=phred${phred} \
				  --log=${name}.umitools.log 


else
umi_tools extract --bc-pattern='${barcode_pattern1}' \
                  --log=${name}.umitools.log \
                  --extract-method=regex \
                  --stdin ${file1} \
                  --stdout result/${name}.fastq.gz \
				  --quality-filter-threshold=${UMIqualityFilterThreshold} \
				  --quality-encoding=phred${phred}
	if [ "${remove_duplicates_based_on_UMI}" == "true" ]; then		  
        mv result/${name}.fastq.gz  result/${name}_umitools.fastq.gz && gunzip result/${name}_umitools.fastq.gz
        ## only checks last part of the underscore splitted header for UMI
        awk '(NR%4==1){name=\$1;header=\$0;len=split(name,umiAr,"_");umi=umiAr[len];} (NR%4==2){total++;if(a[umi]!=1){nondup++;a[umi]=1;  print header;print;getline; print; getline; print;}} END{print FILENAME"\\t"total"\\t"nondup > "${name}.dedup.log"}' result/${name}_umitools.fastq > result/${name}.fastq
        rm result/${name}_umitools.fastq
        gzip result/${name}.fastq
	fi			  
fi
"""

}
}


//* params.run_Trimmer =   "no"   //* @dropdown @options:"yes","no" @show_settings:"Trimmer"
//* @style @multicolumn:{trim_length_5prime,trim_length_3prime}, {trim_length_5prime_R1,trim_length_3prime_R1}, {trim_length_5prime_R2,trim_length_3prime_R2} @condition:{single_or_paired_end_reads="single", trim_length_5prime,trim_length_3prime}, {single_or_paired_end_reads="pair", trim_length_5prime_R1,trim_length_3prime_R1,trim_length_5prime_R2,trim_length_3prime_R2}

//* autofill
//* platform
//* platform
//* autofill
if (!((params.run_Trimmer && (params.run_Trimmer == "yes")) || !params.run_Trimmer)){
g186_23_reads00_g186_19.set{g186_19_reads00_g186_20}
g186_19_log_file10_g186_21 = Channel.empty()
} else {


process Adapter_Trimmer_Quality_Module_Trimmer {

input:
 set val(name), file(reads) from g186_23_reads00_g186_19
 val mate from g_8_1_g186_19

output:
 set val(name), file("reads/*q.gz")  into g186_19_reads00_g186_20
 file "*.log" optional true  into g186_19_log_file10_g186_21

errorStrategy 'retry'

when:
(params.run_Trimmer && (params.run_Trimmer == "yes")) || !params.run_Trimmer

shell:
phred = params.Adapter_Trimmer_Quality_Module_Trimmer.phred
single_or_paired_end_reads = params.Adapter_Trimmer_Quality_Module_Trimmer.single_or_paired_end_reads
trim_length_5prime = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime
trim_length_3prime = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime
trim_length_5prime_R1 = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime_R1
trim_length_3prime_R1 = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime_R1
trim_length_5prime_R2 = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime_R2
trim_length_3prime_R2 = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime_R2
remove_previous_reads = params.Adapter_Trimmer_Quality_Module_Trimmer.remove_previous_reads

nameAll = reads.toString()
nameArray = nameAll.split(' ')
file2 = ""
if (nameAll.contains('.gz')) {
    file1 =  nameArray[0] 
    if (mate == "pair") {file2 =  nameArray[1] }
} 
'''
#!/usr/bin/env perl
 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage; 
 use Cwd qw();
 
runCmd("mkdir reads");
my $file1 = "";
my $file2 = "";
if ("!{mate}" eq "pair") {
    $file1 = "!{file1}";
    $file2 = "!{file2}";
    my $trim1 = "!{trim_length_5prime_R1}:!{trim_length_3prime_R1}";
    my $trim2 = "!{trim_length_5prime_R2}:!{trim_length_3prime_R2}";
    my $len=getLength($file1);
    print "length of $file1: $len\\n";
    trimFiles($file1, $trim1, $len);
    my $len=getLength($file2);
    print "INFO: length of $file2: $len\\n";
    trimFiles($file2, $trim2, $len);
} else {
    $file1 = "!{file1}";
    my $trim1 = "!{trim_length_5prime}:!{trim_length_3prime}";
    my $len=getLength($file1);
    print "INFO: length of file1: $len\\n";
    trimFiles($file1, $trim1, $len);
}
if ("!{remove_previous_reads}" eq "true") {
    my $currpath = Cwd::cwd();
    my @paths = (split '/', $currpath);
    splice(@paths, -2);
    my $workdir= join '/', @paths;
    splice(@paths, -1);
    my $inputsdir= join '/', @paths;
    $inputsdir .= "/inputs";
    print "INFO: inputs reads will be removed if they are located in the workdir inputsdir\\n";
    my @listOfFiles = `readlink -e !{file1} !{file2}`;
    foreach my $targetFile (@listOfFiles){
        if (index($targetFile, $workdir) != -1 || index($targetFile, $inputsdir) != -1) {
            runCmd("rm -f $targetFile");
            print "INFO: $targetFile deleted.\\n";
        }
    }
}



sub trimFiles
{
  my ($file, $trim, $len)=@_;
    my @nts=split(/[,:\\s\\t]+/,$trim);
    my $inpfile="";
    my $com="";
    my $i=1;
    my $outfile="";
    my $param="";
    my $quality="-Q!{phred}";

    if (scalar(@nts)==2)
    {
      $param = "-f ".($nts[0]+1) if (exists($nts[0]) && $nts[0] >= 0 );
      $param .= " -l ".($len-$nts[1]) if (exists($nts[0]) && $nts[1] > 0 );
      $outfile="reads/$file";  
      $com="gunzip -c $file | fastx_trimmer $quality -v $param -z -o $outfile  > !{name}.fastx_trimmer.log" if ((exists($nts[0]) && $nts[0] > 0) || (exists($nts[0]) && $nts[1] > 0 ));
      print "INFO: $com\\n";
      if ($com eq ""){
          print "INFO: Trimmer skipped for $file \\n";
          runCmd("mv $file reads/.");
      } else {
          runCmd("$com");
          print "INFO: Trimmer executed for $file \\n";
      }
    }

    
}


sub getLength
{
   my ($filename)=@_;
   open (IN, $filename);
   my $j=1;
   my $len=0;
   while(my $line=<IN>)
   {
     chomp($line);
     if ($j >50) { last;}
     if ($j%4==0)
     {
        $len=length($line);
     }
     $j++;
   }
   close(IN);
   return $len;
}

sub runCmd {
    my ($com) = @_;
    if ($com eq ""){
		return "";
    }
    my $error = system(@_);
    if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}

'''

}
}



process Adapter_Trimmer_Quality_Module_Trimmer_Removal_Summary {

input:
 file logfile from g186_19_log_file10_g186_21.collect()
 val mate from g_8_1_g186_21

output:
 file "trimmer_summary.tsv"  into g186_21_outputFileTSV00

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %tsvDetail;
my %headerHash;
my %headerText;
my %headerTextDetail;

my $i = 0;
chomp( my $contents = `ls *.log` );

my @files = split( /[\\n]+/, $contents );
foreach my $file (@files) {
    $i++;
    my $mapOrder = "1";
    if ($file =~ /(.*)\\.fastx_trimmer\\.log/){
        $file =~ /(.*)\\.fastx_trimmer\\.log/;
        my $mapper   = "fastx_trimmer";
        my $name = $1;    ##sample name
        push( @header, $mapper );
        my $in;
        my $out;
        chomp( $in =`cat $file | grep 'Input:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $out =`cat $file | grep 'Output:' | awk '{sum+=\\$2} END {print sum}'` );

        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Trimmer" ];
    }
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary          = "trimmer_summary.tsv";
writeFile( $summary,          \\%headerText,       \\%tsv );

sub writeFile {
    my $summary    = $_[0];
    my %headerText = %{ $_[1] };
    my %tsv        = %{ $_[2] };
    open( OUT, ">$summary" );
    print OUT "Sample\\t";
    my @headArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push( @headArr, @{ $headerText{$mapOrder} } );
    }
    my $headArrAll = join( "\\t", @headArr );
    print OUT "$headArrAll\\n";

    foreach my $name ( keys %tsv ) {
        my @rowArr = ();
        for my $mapOrder (@sortedOrderArray) {
            push( @rowArr, @{ $tsv{$name}{ $headerHash{$mapOrder} } } );
        }
        my $rowArrAll = join( "\\t", @rowArr );
        print OUT "$name\\t$rowArrAll\\n";
    }
    close(OUT);
}

'''
}

//* params.run_Quality_Filtering =   "no"   //* @dropdown @options:"yes","no" @show_settings:"Quality_Filtering"
//* @style @multicolumn:{window_size,required_quality}, {leading,trailing,minlen}, {minQuality,minPercent} @condition:{tool="trimmomatic", minlen, trailing, leading, required_quality_for_window_trimming, window_size}, {tool="fastx", minQuality, minPercent}

//* autofill
//* platform
//* platform
//* autofill
if (!((params.run_Quality_Filtering && (params.run_Quality_Filtering == "yes")) || !params.run_Quality_Filtering)){
g186_19_reads00_g186_20.set{g186_20_reads01_g116_31}
g186_20_log_file10_g186_16 = Channel.empty()
} else {


process Adapter_Trimmer_Quality_Module_Quality_Filtering {

input:
 set val(name), file(reads) from g186_19_reads00_g186_20
 val mate from g_8_1_g186_20

output:
 set val(name), file("reads/*.gz")  into g186_20_reads01_g116_31
 file "*.{fastx,trimmomatic}_quality.log" optional true  into g186_20_log_file10_g186_16

when:
(params.run_Quality_Filtering && (params.run_Quality_Filtering == "yes")) || !params.run_Quality_Filtering    

shell:
tool = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.tool
phred = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.phred
window_size = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.window_size
required_quality_for_window_trimming = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.required_quality_for_window_trimming
leading = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.leading
trailing = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.trailing
minlen = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.minlen


// fastx parameters
minQuality = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.minQuality
minPercent = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.minPercent

remove_previous_reads = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.remove_previous_reads

nameAll = reads.toString()
nameArray = nameAll.split(' ')
file2 ="";
if (nameAll.contains('.gz')) {
    file1 =  nameArray[0] 
    if (mate == "pair") {file2 =  nameArray[1]}
} 
'''
#!/usr/bin/env perl
 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage; 
 use Cwd qw();
 
runCmd("mkdir reads unpaired");
my $param = "SLIDINGWINDOW:"."!{window_size}".":"."!{required_quality_for_window_trimming}";
$param.=" LEADING:"."!{leading}";
$param.=" TRAILING:"."!{trailing}";
$param.=" MINLEN:"."!{minlen}";

my $quality="!{phred}";

print "INFO: fastq quality: $quality\\n";
     
if ("!{tool}" eq "trimmomatic") {
    if ("!{mate}" eq "pair") {
        runCmd("trimmomatic PE -phred${quality} !{file1} !{file2} reads/!{name}.1.fastq.gz unpaired/!{name}.1.fastq.unpaired.gz reads/!{name}.2.fastq.gz unpaired/!{name}.1.fastq.unpaired.gz $param 2> !{name}.trimmomatic_quality.log");
    } else {
        runCmd("trimmomatic SE -phred${quality} !{file1} reads/!{name}.fastq.gz $param 2> !{name}.trimmomatic_quality.log");
    }
} elsif ("!{tool}" eq "fastx") {
    if ("!{mate}" eq "pair") {
        print("WARNING: Fastx option is not suitable for paired reads. This step will be skipped.");
        runCmd("mv !{file1} !{file2} reads/.");
    } else {
        runCmd("fastq_quality_filter  -Q $quality -q !{minQuality} -p !{minPercent} -v -i !{file1} -o reads/!{name}.fastq.gz > !{name}.fastx_quality.log");
    }
}
if ("!{remove_previous_reads}" eq "true") {
    my $currpath = Cwd::cwd();
    my @paths = (split '/', $currpath);
    splice(@paths, -2);
    my $workdir= join '/', @paths;
    splice(@paths, -1);
    my $inputsdir= join '/', @paths;
    $inputsdir .= "/inputs";
    print "INFO: inputs reads will be removed if they are located in the workdir inputsdir\\n";
    my @listOfFiles = `readlink -e !{file1} !{file2}`;
    foreach my $targetFile (@listOfFiles){
        if (index($targetFile, $workdir) != -1 || index($targetFile, $inputsdir) != -1) {
            runCmd("rm -f $targetFile");
            print "INFO: $targetFile deleted.\\n";
        }
    }
}

##Subroutines
sub runCmd {
    my ($com) = @_;
    if ($com eq ""){
		return "";
    }
    my $error = system(@_);
    if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}


'''

}
}



process Adapter_Trimmer_Quality_Module_Quality_Filtering_Summary {

input:
 file logfile from g186_20_log_file10_g186_16.collect()
 val mate from g_8_1_g186_16

output:
 file "quality_filter_summary.tsv"  into g186_16_outputFileTSV00

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %headerHash;
my %headerText;

my $i = 0;
chomp( my $contents = `ls *.log` );
my @files = split( /[\\n]+/, $contents );
foreach my $file (@files) {
    $i++;
    my $mapper   = "";
    my $mapOrder = "1";
    if ($file =~ /(.*)\\.fastx_quality\\.log/){
        $mapper   = "fastx";
        $file =~ /(.*)\\.fastx_quality\\.log/;
        my $name = $1;    ##sample name
        push( @header, $mapper );
        my $in;
        my $out;
        chomp( $in =`cat $file | grep 'Input:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $out =`cat $file | grep 'Output:' | awk '{sum+=\\$2} END {print sum}'` );
        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Quality Filtering" ];
    } elsif ($file =~ /(.*)\\.trimmomatic_quality\\.log/){
        $mapper   = "trimmomatic";
        $file =~ /(.*)\\.trimmomatic_quality\\.log/;
        my $name = $1;    ##sample name
        push( @header, $mapper );
        my $in;
        my $out;
        if ( "!{mate}" eq "pair"){
            chomp( $in =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$4} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$7} END {print sum}'` );
        } else {
            chomp( $in =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$3} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$5} END {print sum}'` );
        }
        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Quality Filtering" ];
    }
    
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary          = "quality_filter_summary.tsv";
writeFile( $summary,          \\%headerText,       \\%tsv );

sub writeFile {
    my $summary    = $_[0];
    my %headerText = %{ $_[1] };
    my %tsv        = %{ $_[2] };
    open( OUT, ">$summary" );
    print OUT "Sample\\t";
    my @headArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push( @headArr, @{ $headerText{$mapOrder} } );
    }
    my $headArrAll = join( "\\t", @headArr );
    print OUT "$headArrAll\\n";

    foreach my $name ( keys %tsv ) {
        my @rowArr = ();
        for my $mapOrder (@sortedOrderArray) {
            push( @rowArr, @{ $tsv{$name}{ $headerHash{$mapOrder} } } );
        }
        my $rowArrAll = join( "\\t", @rowArr );
        print OUT "$name\\t$rowArrAll\\n";
    }
    close(OUT);
}

'''
}

//* params.star_index =  ""  //* @input

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 10
    $MEMORY = 50
}
//* platform
//* platform
//* autofill

process STAR_Module_Map_STAR {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}SJ.out.tab$/) "star/$filename"}
input:
 val mate from g_8_0_g116_31
 set val(name), file(reads) from g186_20_reads01_g116_31
 file star_index from g116_26_starIndex02_g116_31

output:
 set val(name), file("${name}Log.final.out")  into g116_31_outputFileOut00_g116_18
 set val(name), file("${name}.flagstat.txt")  into g116_31_outputFileTxt11
 file "${name}Log.out"  into g116_31_logOut22
 set val(name), file("${name}.bam")  into g116_31_mapped_reads30_g116_30
 set val(name), file("${name}SJ.out.tab")  into g116_31_outputFileTab44
 file "${name}Log.progress.out"  into g116_31_progressOut55
 set val(name), file("${name}Aligned.toTranscriptome.out.bam") optional true  into g116_31_transcriptome_bam60_g116_15
 file "${name}Log.final.out"  into g116_31_logFinalOut77

when:
(params.run_STAR && (params.run_STAR == "yes")) || !params.run_STAR

script:
params_STAR = params.STAR_Module_Map_STAR.params_STAR
transcriptomeSAM = ""
if (params.run_Salmon_after_STAR && params.run_Salmon_after_STAR == "yes" && params_STAR.indexOf("--quantMode") < 0){
	transcriptomeSAM = " --quantMode TranscriptomeSAM "
}

"""
STAR --runThreadN ${task.cpus} ${params_STAR} ${transcriptomeSAM} --genomeDir ${star_index} --readFilesCommand zcat --readFilesIn $reads --outFileNamePrefix ${name}
echo "Alignment completed."
if [ ! -e "${name}Aligned.toTranscriptome.out.bam" -a -e "${name}Aligned.toTranscriptome.out.sam" ] ; then
    samtools view -S -b ${name}Aligned.toTranscriptome.out.sam > ${name}Aligned.toTranscriptome.out.bam
elif [ ! -e "${name}Aligned.out.bam" -a -e "${name}Aligned.out.sam" ] ; then
    samtools view -S -b ${name}Aligned.out.sam > ${name}Aligned.out.bam
fi
rm -rf *.sam
if [ -e "${name}Aligned.sortedByCoord.out.bam" ] ; then
    mv ${name}Aligned.sortedByCoord.out.bam ${name}.bam
elif [ -e "${name}Aligned.out.bam" ] ; then
    mv ${name}Aligned.out.bam ${name}.bam
fi

samtools flagstat ${name}.bam > ${name}.flagstat.txt
"""


}


process STAR_Module_STAR_Summary {

input:
 set val(name), file(alignSum) from g116_31_outputFileOut00_g116_18.groupTuple()

output:
 file "*.tsv"  into g116_18_outputFile00_g116_11
 val "star_alignment_sum"  into g116_18_name11_g116_11

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage; 
use Data::Dumper;

my %tsv;
my @headers = ();
my $name = '!{name}';


alteredAligned();

my @keys = keys %tsv;
my $summary = "$name"."_star_sum.tsv";
my $header_string = join("\\t", @headers);
`echo "$header_string" > $summary`;
foreach my $key (@keys){
	my $values = join("\\t", @{ $tsv{$key} });
	`echo "$values" >> $summary`;
}


sub alteredAligned
{
	my @files = qw(!{alignSum});
	my $multimappedSum;
	my $alignedSum;
	my $inputCountSum;
	push(@headers, "Sample");
    push(@headers, "Total Reads");
	push(@headers, "Multimapped Reads Aligned (STAR)");
	push(@headers, "Unique Reads Aligned (STAR)");
	foreach my $file (@files){
		my $multimapped;
		my $aligned;
		my $inputCount;
		chomp($inputCount = `cat $file | grep 'Number of input reads' | awk '{sum+=\\$6} END {print sum}'`);
		chomp($aligned = `cat $file | grep 'Uniquely mapped reads number' | awk '{sum+=\\$6} END {print sum}'`);
		chomp($multimapped = `cat $file | grep 'Number of reads mapped to multiple loci' | awk '{sum+=\\$9} END {print sum}'`);
		$multimappedSum += int($multimapped);
        $alignedSum += int($aligned);
        $inputCountSum += int($inputCount);
	}
	$tsv{$name} = [$name, $inputCountSum];
	push(@{$tsv{$name}}, $multimappedSum);
	push(@{$tsv{$name}}, $alignedSum);
}

sub runCommand {
    my ($com) = @_;
    if ($com eq ""){
		return "";
    }
    my $error = system(@_);
    if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}
'''

}


process STAR_Module_merge_tsv_files_with_same_header {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}.tsv$/) "star_summary/$filename"}
input:
 file tsv from g116_18_outputFile00_g116_11.collect()
 val outputFileName from g116_18_name11_g116_11.collect()

output:
 file "${name}.tsv"  into g116_11_outputFileTSV00_g_118

errorStrategy 'retry'
maxRetries 3

script:
name = outputFileName[0]
"""    
awk 'FNR==1 && NR!=1 {  getline; } 1 {print} ' *.tsv > ${name}.tsv
"""
}

g116_11_outputFileTSV00_g_118= g116_11_outputFileTSV00_g_118.ifEmpty([""]) 

//* autofill
//* platform
//* platform
//* autofill

process Alignment_Summary {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /summary_data.tsv$/) "summary/$filename"}
input:
 file starSum from g116_11_outputFileTSV00_g_118

output:
 file "summary_data.tsv"  into g_118_outputFileTSV00

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_rows;
my @seen_cols;
my $ID_header;

chomp(my $contents = `ls *.tsv`);
my @files = split(/[\\n]+/, $contents);
# if sequential_mapping_short_sum.tsv is exist push to the beginning of the array
my $checkSeqMap = "false";
my $checkSeqMapVal = "";
for my $index (reverse 0..$#files) {
    if ( $files[$index] =~ /sequential_mapping/ ) {
        $checkSeqMap = "true";
        $checkSeqMapVal = $files[$index];
        splice(@files, $index, 1, ());
    }
}
if ($checkSeqMap == "true"){
    unshift @files, $checkSeqMapVal;
}
##Merge each file according to array order

foreach my $file (@files){
        open IN,"$file";
        my $line1 = <IN>;
        chomp($line1);
        ( $ID_header, my @header) = ( split("\\t", $line1) );
        push @seen_cols, @header;

        while (my $line=<IN>) {
        chomp($line);
        my ( $ID, @fields ) = ( split("\\t", $line) ); 
        my %this_row;
        @this_row{@header} = @fields;

        #print Dumper \\%this_row;

        foreach my $column (@header) {
            if (! exists $all_rows{$ID}{$column}) {
                $all_rows{$ID}{$column} = $this_row{$column}; 
            }
        }   
    }
    close IN;
}

#print for debugging
#print Dumper \\%all_rows;
#print Dumper \\%seen_cols;

#grab list of column headings we've seen, and order them. 
my @cols_to_print = uniq(@seen_cols);
my $summary = "summary_data.tsv";
open OUT, ">$summary";
print OUT join ("\\t", $ID_header,@cols_to_print),"\\n";
foreach my $key ( keys %all_rows ) { 
    #map iterates all the columns, and gives the value or an empty string. if it's undefined. (prevents errors)
    print OUT join ("\\t", $key, (map { $all_rows{$key}{$_} // '' } @cols_to_print)),"\\n";
}
close OUT;

sub uniq {
    my %seen;
    grep ! $seen{$_}++, @_;
}

'''


}


//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 10
}
//* platform
//* platform
//* autofill

process STAR_Module_Merge_Bam_and_create_sense_antisense {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_sorted.bam.bai$/) "star/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_sorted.bam$/) "star/$filename"}
input:
 set val(oldname), file(bamfiles) from g116_31_mapped_reads30_g116_30.groupTuple()
 val mate from g_8_1_g116_30

output:
 file "*_sorted.bam.bai"  into g116_30_bam_index00
 set val(oldname), file("*_sorted.bam")  into g116_30_bamFile11_g_161, g116_30_bamFile11_g_175, g116_30_bamFile11_g114_131, g116_30_bamFile10_g114_121, g116_30_bamFile10_g114_143

errorStrategy 'retry'
maxRetries 2

shell:
'''
num=$(echo "!{bamfiles.join(" ")}" | awk -F" " '{print NF-1}')
if [ "${num}" -gt 0 ]; then
    samtools merge !{oldname}.bam !{bamfiles.join(" ")} && samtools sort -o !{oldname}_sorted.bam !{oldname}.bam && samtools index !{oldname}_sorted.bam
else
    mv !{bamfiles.join(" ")} !{oldname}.bam 2>/dev/null || true
    samtools sort  -o !{oldname}_sorted.bam !{oldname}.bam && samtools index !{oldname}_sorted.bam
fi


'''
}


process BAM_Analysis_Module_bam_sort_index {

input:
 set val(name), file(bam) from g116_30_bamFile10_g114_143

output:
 set val(name), file("bam/*.bam"), file("bam/*.bam.bai")  into g114_143_bam_bai00_g114_142, g114_143_bam_bai00_g114_134

when:
params.run_BigWig_Conversion == "yes" || params.run_RSeQC == "yes"

script:
nameAll = bam.toString()
if (nameAll.contains('_sorted.bam')) {
    runSamtools = "samtools index ${nameAll}"
    nameFinal = nameAll
} else {
    runSamtools = "samtools sort -o ${name}_sorted.bam $bam && samtools index ${name}_sorted.bam "
    nameFinal = "${name}_sorted.bam"
}

"""
$runSamtools
mkdir -p bam
mv ${name}_sorted.bam ${name}_sorted.bam.bai bam/.
"""
}

//* params.bed =  ""  //* @input
//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 10
}
//* platform
//* platform
//* autofill

process BAM_Analysis_Module_RSeQC {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*$/) "rseqc/$filename"}
input:
 set val(name), file(bam), file(bai) from g114_143_bam_bai00_g114_134
 file bed from g184_54_bed31_g114_134
 val mate from g_8_2_g114_134

output:
 file "*"  into g114_134_outputFileOut08_g_183

container 'quay.io/viascientific/rseqc:1.0'

when:
(params.run_RSeQC && (params.run_RSeQC == "yes")) || !params.run_RSeQC

script:
run_bam_stat = params.BAM_Analysis_Module_RSeQC.run_bam_stat
run_read_distribution = params.BAM_Analysis_Module_RSeQC.run_read_distribution
run_inner_distance = params.BAM_Analysis_Module_RSeQC.run_inner_distance
run_junction_annotation = params.BAM_Analysis_Module_RSeQC.run_junction_annotation
run_junction_saturation = params.BAM_Analysis_Module_RSeQC.run_junction_saturation
//run_geneBody_coverage and run_infer_experiment needs subsampling
run_geneBody_coverage = params.BAM_Analysis_Module_RSeQC.run_geneBody_coverage
run_infer_experiment = params.BAM_Analysis_Module_RSeQC.run_infer_experiment
"""
if [ "$run_bam_stat" == "true" ]; then bam_stat.py  -i ${bam} > ${name}.bam_stat.txt; fi
if [ "$run_read_distribution" == "true" ]; then read_distribution.py  -i ${bam} -r ${bed}> ${name}.read_distribution.out; fi


if [ "$run_infer_experiment" == "true" -o "$run_geneBody_coverage" == "true" ]; then
	numAlignedReads=\$(samtools view -c -F 4 $bam)

	if [ "\$numAlignedReads" -gt 1000000 ]; then
    	echo "Read number is greater than 1000000. Subsampling..."
    	finalRead=1000000
    	fraction=\$(samtools idxstats  $bam | cut -f3 | awk -v ct=\$finalRead 'BEGIN {total=0} {total += \$1} END {print ct/total}')
    	samtools view -b -s \${fraction} $bam > ${name}_sampled.bam
    	samtools index ${name}_sampled.bam
    	if [ "$run_infer_experiment" == "true" ]; then infer_experiment.py -i ${name}_sampled.bam  -r $bed > ${name}; fi
		if [ "$run_geneBody_coverage" == "true" ]; then geneBody_coverage.py -i ${name}_sampled.bam  -r $bed -o ${name}; fi
	else
		if [ "$run_infer_experiment" == "true" ]; then infer_experiment.py -i $bam  -r $bed > ${name}; fi
		if [ "$run_geneBody_coverage" == "true" ]; then geneBody_coverage.py -i $bam  -r $bed -o ${name}; fi
	fi

fi


if [ "${mate}" == "pair" ]; then
	if [ "$run_inner_distance" == "true" ]; then inner_distance.py -i $bam  -r $bed -o ${name}.inner_distance > stdout.txt; fi
	if [ "$run_inner_distance" == "true" ]; then head -n 2 stdout.txt > ${name}.inner_distance_mean.txt; fi
fi
if [ "$run_junction_annotation" == "true" ]; then junction_annotation.py -i $bam  -r $bed -o ${name}.junction_annotation 2> ${name}.junction_annotation.log; fi
if [ "$run_junction_saturation" == "true" ]; then junction_saturation.py -i $bam  -r $bed -o ${name}; fi
if [ -e class.log ] ; then mv class.log ${name}_class.log; fi
if [ -e log.txt ] ; then mv log.txt ${name}_log.txt; fi
if [ -e stdout.txt ] ; then mv stdout.txt ${name}_stdout.txt; fi


"""

}



//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 30
}
//* platform
//* platform
//* autofill

process BAM_Analysis_Module_UCSC_BAM2BigWig_converter {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.bw$/) "bigwig/$filename"}
input:
 set val(name), file(bam), file(bai) from g114_143_bam_bai00_g114_142
 file genomeSizes from g184_54_genomeSizes21_g114_142

output:
 file "*.bw" optional true  into g114_142_outputFileBw00
 file "publish/*.bw" optional true  into g114_142_publishBw11

container 'quay.io/biocontainers/deeptools:3.5.4--pyhdfd78af_1'

when:
params.run_BigWig_Conversion == "yes"

script:
nameAll = bam.toString()
if (nameAll.contains('_sorted.bam')) {
    runSamtools = ""
    nameFinal = nameAll
} else {
    runSamtools = "mv $bam ${name}_sorted.bam "
    nameFinal = "${name}_sorted.bam"
}
deeptools_parameters = params.BAM_Analysis_Module_UCSC_BAM2BigWig_converter.deeptools_parameters


"""
$runSamtools
bamCoverage ${deeptools_parameters}  -b ${nameFinal} -o ${name}.bw 


"""

}

//* params.pdfbox_path =  ""  //* @input
//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 32
}
//* platform
//* platform
//* autofill

process BAM_Analysis_Module_Picard {

input:
 set val(name), file(bam) from g116_30_bamFile10_g114_121

output:
 file "*_metrics"  into g114_121_outputFileOut00_g114_82
 file "results/*.pdf"  into g114_121_outputFilePdf12_g114_82

container 'public.ecr.aws/t4w5x8f2/viascientific/picard:1.0'

when:
(params.run_Picard_CollectMultipleMetrics && (params.run_Picard_CollectMultipleMetrics == "yes")) || !params.run_Picard_CollectMultipleMetrics

script:
"""
picard CollectMultipleMetrics OUTPUT=${name}_multiple.out VALIDATION_STRINGENCY=LENIENT INPUT=${bam}
mkdir results && java -jar ${params.pdfbox_path} PDFMerger *.pdf results/${name}_multi_metrics.pdf
"""
}


process BAM_Analysis_Module_Picard_Summary {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tsv$/) "picard_summary/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tsv$/) "rseqc_summary/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /results\/.*.pdf$/) "picard/$filename"}
input:
 file picardOut from g114_121_outputFileOut00_g114_82.collect()
 val mate from g_8_1_g114_82
 file picardPdf from g114_121_outputFilePdf12_g114_82.collect()

output:
 file "*.tsv"  into g114_82_outputFileTSV00
 file "results/*.pdf"  into g114_82_outputFilePdf11

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage; 
use Data::Dumper;

runCommand("mkdir results && mv *.pdf results/. ");

my $indir = $ENV{'PWD'};
my $outd = $ENV{'PWD'};
my @files = ();
my @outtypes = ("CollectRnaSeqMetrics", "alignment_summary_metrics", "base_distribution_by_cycle_metrics", "insert_size_metrics", "quality_by_cycle_metrics", "quality_distribution_metrics" );

foreach my $outtype (@outtypes)
{
my $ext="_multiple.out";
$ext.=".$outtype" if ($outtype ne "CollectRnaSeqMetrics");
@files = <$indir/*$ext>;

my @rowheaders=();
my @libs=();
my %metricvals=();
my %histvals=();

my $pdffile="";
my $libname="";
foreach my $d (@files){
  my $libname=basename($d, $ext);
  print $libname."\\n";
  push(@libs, $libname); 
  getMetricVals($d, $libname, \\%metricvals, \\%histvals, \\@rowheaders);
}

my $sizemetrics = keys %metricvals;
write_results("$outd/$outtype.stats.tsv", \\@libs,\\%metricvals, \\@rowheaders, "metric") if ($sizemetrics>0);
my $sizehist = keys %histvals;
write_results("$outd/$outtype.hist.tsv", \\@libs,\\%histvals, "none", "nt") if ($sizehist>0);

}

sub write_results
{
  my ($outfile, $libs, $vals, $rowheaders, $name )=@_;
  open(OUT, ">$outfile");
  print OUT "$name\\t".join("\\t", @{$libs})."\\n";
  my $size=0;
  $size=scalar(@{${$vals}{${$libs}[0]}}) if(exists ${$libs}[0] and exists ${$vals}{${$libs}[0]} );
  
  for (my $i=0; $i<$size;$i++)
  { 
    my $rowname=$i;
    $rowname = ${$rowheaders}[$i] if ($name=~/metric/);
    print OUT $rowname;
    foreach my $lib (@{$libs})
    {
      print OUT "\\t".${${$vals}{$lib}}[$i];
    } 
    print OUT "\\n";
  }
  close(OUT);
}

sub getMetricVals{
  my ($filename, $libname, $metricvals, $histvals,$rowheaders)=@_;
  if (-e $filename){
     my $nextisheader=0;
     my $nextisvals=0;
     my $nexthist=0;
     open(IN, $filename);
     while(my $line=<IN>)
     {
       chomp($line);
       @{$rowheaders}=split(/\\t/, $line) if ($nextisheader && !scalar(@{$rowheaders})); 
       if ($nextisvals) {
         @{${$metricvals}{$libname}}=split(/\\t/, $line);
         $nextisvals=0;
       }
       if($nexthist){
          my @vals=split(/[\\s\\t]+/,$line); 
          push(@{${$histvals}{$libname}}, $vals[1]) if (exists $vals[1]);
       }
       $nextisvals=1 if ($nextisheader); $nextisheader=0;
       $nextisheader=1 if ($line=~/METRICS CLASS/);
       $nexthist=1 if ($line=~/normalized_position/);
     } 
  }
  
}


sub runCommand {
	my ($com) = @_;
	if ($com eq ""){
		return "";
    }
    my $error = system(@_);
	if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}
'''

}

igv_extention_factor = params.BAM_Analysis_Module_IGV_BAM2TDF_converter.igv_extention_factor
igv_window_size = params.BAM_Analysis_Module_IGV_BAM2TDF_converter.igv_window_size

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 24
}
//* platform
//* platform
//* autofill

process BAM_Analysis_Module_IGV_BAM2TDF_converter {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tdf$/) "igvtools_result/$filename"}
input:
 val mate from g_8_0_g114_131
 set val(name), file(bam) from g116_30_bamFile11_g114_131
 file genomeSizes from g184_54_genomeSizes22_g114_131

output:
 file "*.tdf"  into g114_131_outputFileOut00

when:
(params.run_IGV_TDF_Conversion && (params.run_IGV_TDF_Conversion == "yes")) || !params.run_IGV_TDF_Conversion

script:
pairedText = (params.nucleicAcidType == "dna" && mate == "pair") ? " --pairs " : ""
nameAll = bam.toString()
if (nameAll.contains('_sorted.bam')) {
    runSamtools = "samtools index ${nameAll}"
    nameFinal = nameAll
} else {
    runSamtools = "samtools sort -o ${name}_sorted.bam $bam && samtools index ${name}_sorted.bam "
    nameFinal = "${name}_sorted.bam"
}
"""
$runSamtools
igvtools count -w ${igv_window_size} -e ${igv_extention_factor} ${pairedText} ${nameFinal} ${name}.tdf ${genomeSizes}
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 52
}
//* autofill

process ESAT_bulk_XL_copy {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*umi.distributions.txt$/) "ESAT_withMultipleMapping_UMI/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.log$/) "ESAT_withMultipleMapping_log/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.gene.txt$/) "ESAT_withMultipleMapping_gene/$filename"}
input:
 file trans2gene_filepath from g_125_0_g_175
 set val(name), file(sorted_bam) from g116_30_bamFile11_g_175

output:
 set val(name), file("*umi.distributions.txt")  into g_175_UMI_distributions00_g_176
 file "*.log"  into g_175_log_file11
 file "*.gene.txt"  into g_175_txtFile22

"""
find  -name "*.bam" | awk '{print "${name}\t"\$1 }' > filelist.txt

java \
-Xmx50g \
-jar  /usr/local/bin/esat.v0.1_09.09.16_24.18.umihack.jar \
-alignments filelist.txt \
-out ${name}_esat.txt \
-geneMapping ${trans2gene_filepath} \
-task score3p \
-wLen 100 \
-wOlap 50 \
-wExt 1000 \
-sigTest .01 \
-multimap normal \
-scPrep \
-umiMin 1
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 30
}
//* platform
//* platform
//* autofill

process UMI_Trim_bulk_XL_copy {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_umi_count.txt$/) "umiCount_withMultipleMapping/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_reads_count.txt$/) "readsCount_withMultipleMapping/$filename"}
input:
 set val(name), file(umi_dist) from g_175_UMI_distributions00_g_176.groupTuple()

output:
 set val(name), file("*_umi_count.txt")  into g_176_UMI_clean00
 set val(name), file("*_reads_count.txt")  into g_176_count_file11

"""
path=\$(which cleanLowEndUmis_withRaw_forBulk.py) && cp \$path .

python cleanLowEndUmis_withRaw_forBulk.py \
-i ${umi_dist} \
-o ${name}_umi_count.txt \
-r ${name}_reads_count.txt \
-u 1
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 52
}
//* autofill

process ESAT_bulk_XL {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*umi.distributions.txt$/) "umi_distribution/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.log$/) "ESAT_log/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.gene.txt$/) "ESAT_gene_count/$filename"}
input:
 file trans2gene_filepath from g_125_0_g_161
 set val(name), file(sorted_bam) from g116_30_bamFile11_g_161

output:
 set val(name), file("*umi.distributions.txt")  into g_161_UMI_distributions00_g_137
 file "*.log"  into g_161_log_file11
 file "*.gene.txt"  into g_161_txtFile22

"""
find  -name "*.bam" | awk '{print "${name}\t"\$1 }' > filelist.txt

java \
-Xmx50g \
-jar  /usr/local/bin/esat.v0.1_09.09.16_24.18.umihack.jar \
-alignments filelist.txt \
-out ${name}_esat.txt \
-geneMapping ${trans2gene_filepath} \
-task score3p \
-wLen 100 \
-wOlap 50 \
-wExt 1000 \
-sigTest .01 \
-multimap proper \
-scPrep \
-umiMin 1
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 30
}
//* platform
//* platform
//* autofill

process UMI_Trim_bulk_XL {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_umi_count.txt$/) "umiCount/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_reads_count.txt$/) "readCount/$filename"}
input:
 set val(name), file(umi_dist) from g_161_UMI_distributions00_g_137.groupTuple()

output:
 set val(name), file("*_umi_count.txt")  into g_137_UMI_clean00
 set val(name), file("*_reads_count.txt")  into g_137_count_file11

"""
path=\$(which cleanLowEndUmis_withRaw_forBulk.py) && cp \$path .

python cleanLowEndUmis_withRaw_forBulk.py \
-i ${umi_dist} \
-o ${name}_umi_count.txt \
-r ${name}_reads_count.txt \
-u 1
"""
}


//* autofill
//* platform
//* platform
//* autofill

process STAR_Module_merge_transcriptome_bam {

input:
 set val(oldname), file(bamfiles) from g116_31_transcriptome_bam60_g116_15.groupTuple()

output:
 set val(oldname), file("${oldname}.bam")  into g116_15_merged_bams00
 set val(oldname), file("*_sorted*bai")  into g116_15_bam_index11
 set val(oldname), file("*_sorted*bam")  into g116_15_sorted_bam22

shell:
'''
num=$(echo "!{bamfiles.join(" ")}" | awk -F" " '{print NF-1}')
if [ "${num}" -gt 0 ]; then
    samtools merge !{oldname}.bam !{bamfiles.join(" ")} && samtools sort -o !{oldname}_sorted.bam !{oldname}.bam && samtools index !{oldname}_sorted.bam
else
    mv !{bamfiles.join(" ")} !{oldname}.bam 2>/dev/null || true
    samtools sort  -o !{oldname}_sorted.bam !{oldname}.bam && samtools index !{oldname}_sorted.bam
fi
'''
}


process Adapter_Trimmer_Quality_Module_Umitools_Summary {

input:
 file logfile from g186_23_log_file10_g186_24.collect()
 val mate from g_8_1_g186_24

output:
 file "umitools_summary.tsv"  into g186_24_outputFileTSV00

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %tsvDetail;
my %headerHash;
my %headerText;
my %headerTextDetail;

my $i = 0;
chomp( my $contents = `ls *.log` );

my @files = split( /[\\n]+/, $contents );
foreach my $file (@files) {
    $i++;
    my $mapOrder = "1";
    if ($file =~ /(.*)\\.umitools\\.log/){
        $file =~ /(.*)\\.umitools\\.log/;
        my $mapper   = "umitools";
        my $name = $1;    ##sample name
        push( @header, $mapper );
        my $in;
        my $out;
        my $dedupout;
        chomp( $in =`cat $file | grep 'INFO Input Reads:' | awk '{sum=\\$6} END {print sum}'` );
        chomp( $out =`cat $file | grep 'INFO Reads output:' | awk '{sum=\\$6} END {print sum}'` );
        my $deduplog = $name.".dedup.log";
        $headerHash{$mapOrder} = $mapper;
        if (-e $deduplog) {
            print "dedup log found\\n";
            chomp( $dedupout =`cat $deduplog | grep '$name' | awk '{sum=\\$3} END {print sum}'` );
            $tsv{$name}{$mapper} = [ $in, $out, $dedupout];
            $headerText{$mapOrder} = [ "Total Reads", "Reads After Umiextract", "Reads After Deduplication" ]; 
        } else {
            $tsv{$name}{$mapper} = [ $in, $out ];
            $headerText{$mapOrder} = [ "Total Reads", "Reads After Umiextract" ]; 
        }
        
        
    }
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary          = "umitools_summary.tsv";
writeFile( $summary,          \\%headerText,       \\%tsv );

sub writeFile {
    my $summary    = $_[0];
    my %headerText = %{ $_[1] };
    my %tsv        = %{ $_[2] };
    open( OUT, ">$summary" );
    print OUT "Sample\\t";
    my @headArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push( @headArr, @{ $headerText{$mapOrder} } );
    }
    my $headArrAll = join( "\\t", @headArr );
    print OUT "$headArrAll\\n";

    foreach my $name ( keys %tsv ) {
        my @rowArr = ();
        for my $mapOrder (@sortedOrderArray) {
            push( @rowArr, @{ $tsv{$name}{ $headerHash{$mapOrder} } } );
        }
        my $rowArrAll = join( "\\t", @rowArr );
        print OUT "$name\\t$rowArrAll\\n";
    }
    close(OUT);
}

'''
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 10
}
//* platform
//* platform
//* autofill

process Adapter_Trimmer_Quality_Module_FastQC_after_Adapter_Removal {

input:
 val mate from g_8_0_g186_31
 set val(name), file(reads) from g186_18_reads01_g186_31

output:
 file '*.{html,zip}'  into g186_31_FastQCout00

when:
(params.run_FastQC && params.run_FastQC == "yes" && params.run_Adapter_Removal && params.run_Adapter_Removal == "yes")

script:
"""
fastqc ${reads} 
"""
}


process Adapter_Trimmer_Quality_Module_Adapter_Removal_Summary {

input:
 file logfile from g186_18_log_file10_g186_11.collect()
 val mate from g_8_1_g186_11

output:
 file "adapter_removal_summary.tsv"  into g186_11_outputFileTSV00
 file "adapter_removal_detailed_summary.tsv" optional true  into g186_11_outputFile11

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %tsvDetail;
my %headerHash;
my %headerText;
my %headerTextDetail;

my $i = 0;
chomp( my $contents = `ls *.log` );

my @files = split( /[\\n]+/, $contents );
foreach my $file (@files) {
    $i++;
    my $mapOrder = "1";
    if ($file =~ /(.*)\\.fastx\\.log/){
        $file =~ /(.*)\\.fastx\\.log/;
        my $mapper   = "fastx";
        my $name = $1;    ##sample name
        push( @header, $mapper );

        my $in;
        my $out;
        my $tooshort;
        my $adapteronly;
        my $noncliped;
        my $Nreads;

        chomp( $in =`cat $file | grep 'Input:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $out =`cat $file | grep 'Output:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $tooshort =`cat $file | grep 'too-short reads' | awk '{sum+=\\$2} END {print sum}'`);
        chomp( $adapteronly =`cat $file | grep 'adapter-only reads' | awk '{sum+=\\$2} END {print sum}'`);
        chomp( $noncliped =`cat $file | grep 'non-clipped reads.' | awk '{sum+=\\$2} END {print sum}'`);
        chomp( $Nreads =`cat $file | grep 'N reads.' | awk '{sum+=\\$2} END {print sum}'` );

        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Adapter Removal" ];
        $tsvDetail{$name}{$mapper} = [ $in, $tooshort, $adapteronly, $noncliped, $Nreads, $out ];
        $headerTextDetail{$mapOrder} = ["Total Reads","Too-short reads","Adapter-only reads","Non-clipped reads","N reads","Reads After Adapter Removal"];
    } elsif ($file =~ /(.*)\\.trimmomatic\\.log/){
        $file =~ /(.*)\\.trimmomatic\\.log/;
        my $mapper   = "trimmomatic";
        my $name = $1;    ##sample name
        push( @header, $mapper );
        
        my $in;
        my $out;

        if ( "!{mate}" eq "pair"){
            chomp( $in =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$4} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$7} END {print sum}'` );
        } else {
            chomp( $in =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$3} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$5} END {print sum}'` );
        }
        


        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Adapter Removal" ];
        
    }
    
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary          = "adapter_removal_summary.tsv";
my $detailed_summary = "adapter_removal_detailed_summary.tsv";
writeFile( $summary,          \\%headerText,       \\%tsv );
if (%headerTextDetail){
    writeFile( $detailed_summary, \\%headerTextDetail, \\%tsvDetail );  
}

sub writeFile {
    my $summary    = $_[0];
    my %headerText = %{ $_[1] };
    my %tsv        = %{ $_[2] };
    open( OUT, ">$summary" );
    print OUT "Sample\\t";
    my @headArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push( @headArr, @{ $headerText{$mapOrder} } );
    }
    my $headArrAll = join( "\\t", @headArr );
    print OUT "$headArrAll\\n";

    foreach my $name ( keys %tsv ) {
        my @rowArr = ();
        for my $mapOrder (@sortedOrderArray) {
            push( @rowArr, @{ $tsv{$name}{ $headerHash{$mapOrder} } } );
        }
        my $rowArrAll = join( "\\t", @rowArr );
        print OUT "$name\\t$rowArrAll\\n";
    }
    close(OUT);
}

'''
}


//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 10
}
//* platform
//* platform
//* autofill

process Adapter_Trimmer_Quality_Module_FastQC {

input:
 val mate from g_8_0_g186_28
 set val(name), file(reads) from g_166_valid_fastq01_g186_28

output:
 file '*.{html,zip}'  into g186_28_FastQCout04_g_183

when:
(params.run_FastQC && (params.run_FastQC == "yes"))

script:
"""
fastqc ${reads} 
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 24
}
//* platform
//* platform
//* autofill

process MultiQC {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /multiqc_report.html$/) "multiQC/$filename"}
input:
 file "fastqc/*" from g186_28_FastQCout04_g_183.flatten().toList()
 file "rseqc_star/*" from g114_134_outputFileOut08_g_183.flatten().toList()

output:
 file "multiqc_report.html" optional true  into g_183_outputHTML00
 file "*" optional true  into g_183_outputDir11


script:
multiqc_parameters = params.MultiQC.multiqc_parameters
plots_flat_numseries = params.MultiQC.plots_flat_numseries

"""
multiqc ${multiqc_parameters}  -d -dd 2 --cl_config "plots_flat_numseries: ${plots_flat_numseries}" .

"""


}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
