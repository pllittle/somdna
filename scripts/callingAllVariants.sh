#!/bin/sh

[ ! -z $srcPL_dnaCallAll ] && [ $srcPL_dnaCallAll -eq 1 ] && return 0
[ -z "$git_dir" ] && git_dir=$(cd $(dirname $BASH_SOURCE)/../..; pwd)

for repo in baSHic; do
	repo_dir=$git_dir/$repo
	check_array $repo baSHic && tmp_url=https://github.com/pllittle/$repo.git
	
	while true; do
		if [ ! -d "$repo_dir" ]; then
			cd "$git_dir"
			git clone "$tmp_url" >&2
			[ $? -eq 0 ] && break
		else
			cd "$repo_dir"
			git pull >&2
			[ $? -eq 0 ] && break
		fi
		echo -e "Some error in cloning $repo, contact pllittle" >&2 && return 1
	done
	
done

for fn in genomic; do
	. $git_dir/somdna/scripts/$fn.sh
	[ $? -eq 0 ] && continue
	echo -e "Error src-ing $fn.sh" >&2 && return 1
done

run_strelka2_soma(){
	local gatk_dir strelka_dir nbam tbam ref out_dir ncores
	local regions confirm config_fn var_dir status
	
	confirm=0
	while [ ! -z $1 ]; do
		case $1 in
			-s | --strelka_dir )
				shift
				strelka_dir="$1"
				;;
			-g | --gatk_dir )
				shift
				gatk_dir="$1"
				;;
			-t | --tbam )
				shift
				tbam="$1"
				;;
			-n | --nbam )
				shift
				nbam="$1"
				;;
			-r | --ref )
				shift
				ref="$1"
				;;
			-o | --out_dir )
				shift
				out_dir="$1"
				;;
			-c | --ncores )
				shift
				ncores="$1"
				;;
			--regions )
				shift
				regions="$1"
				;;
			--confirm )
				confirm=1
				;;
			--config_fn )
				shift
				config_fn="$1"
				;;
		esac
		shift
	done
	
	# Check inputs
	[ -z "$strelka_dir" ] && echo "Add -s <Strelka2 dir>" >&2 && return 1
	[ -z "$gatk_dir" ] 		&& echo "Add -g <GATK dir>" >&2 && return 1
	[ -z "$tbam" ] 				&& echo "Add -t <tumor bam>" >&2 && return 1
	[ -z "$nbam" ] 				&& echo "Add -n <normal bam>" >&2 && return 1
	[ -z "$ref" ] 				&& echo "Add -r <reference genome>" >&2 && return 1
	[ -z "$out_dir" ] 		&& echo "Add -o <output dir>" >&2 && return 1
	[ -z "$ncores" ] 			&& echo "Add -c <number of threads/cores>" >&2 && return 1
	
	[ ! $($strelka_dir/bin/configureStrelkaSomaticWorkflow.py -h > /dev/null; echo $?) -eq 0 ] \
		&& echo "Error: Strelka2 not properly installed or environment isn't setup yet" >&2 \
		&& return 1
	
	new_mkdir "$out_dir"
	var_dir="$out_dir/results/variants"
	[ -s "$out_dir/somatic.vcf.gz" ] && return 0
	
	# Prepare regions
	[ ! -z "$regions" ] \
		&& regions=$(echo $regions | sed 's|^|,|g' | sed 's|,| --region |g' | sed 's|^ ||g')
	
	# Clear out any older stuff
	new_rm $out_dir/results $out_dir/workspace \
		$out_dir/runWorkflow* $out_dir/workflow*
	
	# Configure strelka stuff
	local cmd
	cmd="$strelka_dir/bin/configureStrelkaSomaticWorkflow.py"
	cmd="$cmd --normalBam $nbam --tumorBam $tbam --referenceFasta $ref"
	cmd="$cmd --disableEVS --exome"
	[ ! -z $config_fn ] && cmd="$cmd --config $config_fn"
	[ ! -z "$regions" ] && cmd="$cmd $regions"
	cmd="$cmd --runDir $out_dir >&2"
	
	if [ $confirm -eq 1 ]; then
		local resp
		echo -e "`date`: Strelka Command:\n\n$cmd\n\n" >&2
		make_menu -y -p "Does this look good?"; read resp
		[ -z $resp ] && return 1
		[ ! -z $resp ] && [ ! $resp -eq 1 ] && return 1
	fi
	
	# Run variant calling
	if [ ! -f $var_dir/somatic.snvs.vcf.gz ] || [ ! -f $var_dir/somatic.indels.vcf.gz ]; then
		eval $cmd
		status=$?
		[ ! $status -eq 0 ] && echo -e "`date`: Strelka configure error" >&2 && return 1
	
		export OMP_NUM_THREADS=$ncores
		$out_dir/runWorkflow.py -m local -j $ncores >&2
		status=$?
		[ $status -eq 0 ] \
			&& echo -e "`date`: Strelka2 completed" >&2 \
			|| echo -e "`date`: Strelka2 incomplete" >&2
		[ ! $status -eq 0 ] && return 1
	fi
	
	# Merge SNV and INDEL vcfs
	export OMP_NUM_THREADS=1
	$gatk_dir/gatk MergeVcfs \
		--INPUT $var_dir/somatic.snvs.vcf.gz \
		--INPUT $var_dir/somatic.indels.vcf.gz \
		--OUTPUT $out_dir/somatic.vcf.gz >&2
	status=$?
	[ ! $status -eq 0 ] && echo -e "`date`: Error with MergeVcfs" >&2 && return 1
	if [ ! -s $out_dir/somatic.vcf.gz ]; then
		new_rm $out_dir/somatic.vcf.gz
		echo -e "`date`: Empty vcf" >&2 && return 1
	fi
	
	# Clean up Strelka stuff
	new_rm $out_dir/results $out_dir/workspace \
		$out_dir/runWorkflow* $out_dir/workflow*
	
	return 0
	
}
run_strelka2_germ(){
	local stk2_dir bams bams2 ref out_dir ncores
	local regions name out_fn
	
	while [ ! -z $1 ]; do
		case $1 in
			-b | --bams )
				shift
				bams="$1" # comma-delimited array of bams
				;;
			-c | --ncores )
				shift
				ncores="$1"
				;;
			-g | --regions )
				shift
				regions="$1"
				;;
			-n | --name )
				shift
				name="$1"
				;;
			-o | --out_dir )
				shift
				out_dir="$1"
				;;
			-r | --ref )
				shift
				ref="$1"
				;;
			-s | --stk2_dir )
				shift
				stk2_dir="$1"
				;;
		esac
		shift
	done
	
	[ -z $stk2_dir ] 	&& echo "Add -s <Strelka2 directory>" >&2 && return 1
	[ -z $bams ] 			&& echo "Add -b <common delimited array of bams>" >&2 && return 1
	[ -z $name ] 			&& echo "Add -n <name for output file>" >&2 && return 1
	[ -z $ref ] 			&& echo "Add -r <fasta reference fn>" >&2 && return 1
	[ -z $out_dir ] 	&& echo "Add -o <out_dir>" >&2 && return 1
	[ -z $ncores ] 		&& echo "Add -c <ncores>" >&2 && return 1
	
	new_mkdir $out_dir
	out_fn=$out_dir/$name.variants.vcf.gz
	[ -f $out_fn ] && [ -f $out_fn.tbi ] \
		&& return 0
	
	# Clear out any older stuff
	new_rm $out_dir/results $out_dir/workspace \
		$out_dir/runWorkflow* $out_dir/workflow*
	
	# Prepare regions
	[ ! -z "$regions" ] \
		&& regions=$(echo $regions | sed 's|^|,|g' | sed 's|,| --region |g' | sed 's|^ ||g')
	
	# Configuration
	local cmd
	bams2=`echo $bams | sed 's|^|,|g' | sed 's|,| --bam |g' | sed 's|^ ||g'`
	cmd="$stk2_dir/bin/configureStrelkaGermlineWorkflow.py"
	cmd="$cmd $bams2 --referenceFasta $ref --disableEVS --exome"
	[ ! -z "$regions" ] && cmd="$cmd $regions"
	cmd="$cmd --runDir $out_dir"
	# echo -e "Strelka Command:\n\n\t$cmd\n\n" >&2
	eval $cmd >&2
	[ ! $? -eq 0 ] && echo "Error in configuration" >&2 && return 1
	
	# execution on a single local machine with ncores parallel jobs
	export OMP_NUM_THREADS=$ncores
	$out_dir/runWorkflow.py -m local -j $ncores >&2
	[ ! $? -eq 0 ] && echo "There was an error" >&2 && return 1
	export OMP_NUM_THREADS=1
	
	# Clean up
	mv $out_dir/results/variants/variants.vcf.gz \
		$out_dir/$name.variants.vcf.gz
	mv $out_dir/results/variants/variants.vcf.gz.tbi \
		$out_dir/$name.variants.vcf.gz.tbi
	new_rm $out_dir/results $out_dir/workspace \
		$out_dir/runWorkflow* $out_dir/workflow.*
	
	return 0
	
}
down_cosmic(){
	local genome version email pw url authstr downurl out_fn out_dir
	local tmp_fn=~/.down
	
	while [ ! -z $1 ]; do
		case $1 in
			-g | --genome )
				shift
				genome="$1"
				;;
			-v | --version )
				shift
				version="$1"
				;;
			-o | --out_dir )
				shift
				out_dir="$1"
				;;
		esac
		shift
	done
	
	[ -z $genome ] 	&& echo "Add -g <genome>, like GRCh37 or GRCh38" >&2 && return 1
	[ -z $version ] && echo "Add -v <version>, like 98 or 99" >&2 && return 1
	[ -z $out_dir ] && echo "Add -o <out dir>" >&2 && return 1
	new_mkdir $out_dir
	cd $out_dir
	
	url=https://cancer.sanger.ac.uk/cosmic/file_download
	url=$url/$genome/cosmic/v$version/VCF/CosmicCodingMuts.vcf.gz
	# out_fn=`echo $url | sed 's|/|\n|g' | tail -n 1`
	out_fn=CosmicCodingMuts_${genome}_v${version}.vcf.gz
	
	if [ ! -f $out_fn ]; then
		make_menu -p "COSMIC/Sanger Email? (e.g. abc@gmail.com)"; read email
		make_menu -p "COSMIC/Sanger Password?"; read -s pw
		echo >&2
		authstr=`echo -e "${email}:${pw}" | base64`
		curl -H "Authorization: Basic ${authstr}" ${url} > $tmp_fn
		downurl=`cat $tmp_fn | tail -n 1 | sed 's|"||g' \
			| cut -d ':' --complement -f1 \
			| sed 's|}$||g'`
		rm $tmp_fn
		echo $downurl > cosmic_downurl.txt
		curl "${downurl}" -o $out_fn >&2
		[ ! $? -eq 0 ] && echo -e "${red}Error in COSMIC download${NC}" >&2 && return 1
		new_rm cosmic_downurl.txt
	fi
	
}
get_COSMIC_canonical(){
	local genome version cosm_dir
	local hts_dir cosmic_fn
	
	while [ ! -z $1 ]; do
		case $1 in
			-g | --genome )
				shift
				genome="$1"
				;;
			-v | --version )
				shift
				version="$1"
				;;
			-c | --cosm_dir )
				shift
				cosm_dir="$1"
				;;
			-h | --hts_dir )
				shift
				hts_dir="$1"
				;;
		esac
		shift
	done
	
	[ -z $genome ] 		&& echo "Add -g <genome>, e.g. GRCh37/GRCh38" >&2 && return 1
	[ -z $version ] 	&& echo "Add -v <version>, e.g. 94/95/99" >&2 && return 1
	[ -z $hts_dir ] 	&& echo "Add -h <hts_dir>" >&2 && return 1
	[ -z $cosm_dir ] 	&& echo "Add -c <COSMIC dir>" >&2 && return 1
	
	# down_cosmic -g $genome -v $version -o $cosm_dir
	cosmic_fn=$cosm_dir/CosmicCodingMuts_${genome}_v${version}
	
	# Check dependencies
	for fn in bgzip tabix htsfile; do
			[ ! -f "$hts_dir/bin/$fn" ] \
				&& echo "Install htslib" >&2 && return 1
			$hts_dir/bin/$fn --help > /dev/null
			[ $? -eq 0 ] && continue
			echo "Load htslib and dependencies" >&2 && return 1
	done
	
	if [ ! -f ${cosmic_fn}_canonical.vcf.gz ] \
		|| [ ! -f ${cosmic_fn}_canonical.vcf.gz.tbi ]; then
		
		[ ! -f $cosmic_fn.vcf.gz ] \
			&& down_cosmic -g $genome -v $version -o $cosm_dir
		
		echo -e "`date`: Removing some rows" >&2
		zgrep -v "GENE=.*_ENST[0-9]*;" $cosmic_fn.vcf.gz \
			> ${cosmic_fn}_canonical.vcf
		
		echo -e "`date`: Running bgzip ..." >&2
		$hts_dir/bin/bgzip -c ${cosmic_fn}_canonical.vcf \
			> ${cosmic_fn}_canonical.vcf.gz
		[ ! $? -eq 0 ] && echo "Error with bgzip" >&2 && return 1
		
		echo -e "`date`: Running tabix ..." >&2
		$hts_dir/bin/tabix -p vcf ${cosmic_fn}_canonical.vcf.gz
		[ ! $? -eq 0 ] && echo "Error with tabix" >&2 && return 1
		new_rm ${cosmic_fn}_canonical.vcf $cosmic_fn.vcf.gz
		
		echo -e "`date`: Finished downloading/processing COSMIC file for VEP" >&2
		
	else
		echo -e "`date`: CosmicCodingMuts_${genome}_v${version} file already available ^_^" >&2
		
	fi
	
	return 0
}
run_VEP(){
	local fasta_fn vep_dir genome status vep_rel cmd
	local input_fn output_fn vep_fields cosmic_fn ncores
	local cache_dir species cache_type type_dir genome_dir
	local gnomad_fn
	
	local vep_cache0 vep_cache_dir
	
	ncores=1
	while [ ! -z "$1" ]; do
		case $1 in
			-a | --gnomad_fn )
				shift
				gnomad_fn="$1"
				;;
			-c | --cosmic_fn )
				shift
				cosmic_fn="$1"
				;;
			-f | --fasta_fn )
				shift
				fasta_fn="$1"
				;;
			-g | --genome )
				shift
				genome="$1"
				;;
			-h | --cache_dir )
				shift
				cache_dir="$1"
				;;
			-i | --input_fn )
				shift
				input_fn="$1"
				;;
			-n | --ncores )
				shift
				ncores="$1"
				;;
			-o | --output_fn )
				shift
				output_fn="$1"
				;;
			-r | --vep_rel )
				shift
				vep_rel="$1"
				;;
			-s | --species )
				shift
				species="$1"
				;;
			-t | --cache_type )
				shift
				cache_type="$1"
				;;
			-v | --vep_dir )
				shift
				vep_dir="$1"
				;;
		esac
		shift
	done
	
	# Check inputs
	[ -z "$cosmic_fn" ] && echo "Add -c <cosmic_fn>" >&2 && return 1
	[ -z "$gnomad_fn" ] && echo "Add -a <gnomad_fn>" >&2 && return 1
	[ -z "$fasta_fn" ] 	&& echo "Add -f <fasta_fn>" >&2 && return 1
	[ -z "$genome" ] 		&& echo "Add -g <genome, e.g. GRCh37>" >&2 && return 1
	[ -z "$input_fn" ] 	&& echo "Add -i <input_fn>" >&2 && return 1
	[ -z "$output_fn" ] && echo "Add -o <output_fn>" >&2 && return 1
	[ -z "$vep_dir" ] 	&& echo "Add -v <vep_dir>" >&2 && return 1
	[ -z "$vep_rel" ] 	&& echo "Add -r <vep release number, e.g 111>" >&2 && return 1
	[ -z "$cache_dir" ] && echo "Add -h <cache_dir>" >&2 && return 1
	[ -z "$species" ] 	&& species=homo_sapiens
	[ -z "$cache_type" ] && echo "Add -t <cache_type, e.g. vep, refseq, merged>" >&2 && return 1
	
	# Check VEP installed
	[ ! -f $vep_dir/vep ] && echo "Error: VEP missing" >&2 && return 1
	[ ! $($vep_dir/vep --help > /dev/null; echo $?) -eq 0 ] \
		&& echo "Error: VEP not installed or environment not setup" >&2 \
		&& return 1
	
	# Check files
	[ ! -f "$gnomad_fn" ] && echo -e "$gnomad_fn missing" >&2 && return 1
	[ ! -f "$cosmic_fn" ] && echo -e "$cosmic_fn missing" >&2 && return 1
	
	# Check cache_dir
	[ ! -d "$cache_dir" ] && echo -e "$cache_dir missing" >&2 && return 1
	
	# Check cache_type and dir
	check_array "$cache_type" vep refseq merged
	[ ! $? -eq 0 ] && echo "cache_type should be vep, refseq, or merged" >&2 && return 1
	[ "$cache_type" == "vep" ] && type_dir="$cache_dir/$species"
	[ "$cache_type" != "vep" ] && type_dir="$cache_dir/${species}_${cache_type}"
	[ ! -d "$type_dir" ] && echo -e "$type_dir missing" >&2 && return 1
	
	# Check type_dir's genome
	genome_dir=$type_dir/${vep_rel}_${genome}
	[ ! -d "$genome_dir" ] && echo -e "$genome_dir missing" >&2 && return 1
	
	# If output file exists, done
	[ -f $output_fn.gz ] && echo -e "$(date): Final VEP output already exists" >&2 && return 0
	
	echo -e "`date`: Start VEP" >&2
	
	# Run VEP
	vep_fields=IMPACT,Consequence,SYMBOL,HGVSc,HGVSp,AF
	vep_fields="$vep_fields,gnomAD_AF"
	vep_fields="$vep_fields,COSMIC,COSMIC_CNT,COSMIC_LEGACY_ID"
	
	# echo "Editing code here" >&2 && return 1
	
	export OMP_NUM_THREADS=$ncores
	cmd="$vep_dir/vep --format vcf --species $species"
	cmd="$cmd -i $input_fn -o $output_fn --fork $ncores"
	cmd="$cmd --cache --dir_cache $cache_dir --cache_version $vep_rel"
	[ "$cache_type" != "vep" ] && cmd="$cmd --$cache_type"
	cmd="$cmd --assembly $genome --fasta $fasta_fn --force_overwrite"
	cmd="$cmd --no_stats --domains --hgvs --af --vcf"
	# --af_gnomad
	cmd="$cmd --custom file=$cosmic_fn,short_name=COSMIC,format=vcf,type=exact,coords=0,fields=CNT%LEGACY_ID"
	cmd="$cmd --custom file=$gnomad_fn,short_name=gnomAD,format=vcf,type=exact,coords=0,fields=AF"
	cmd="$cmd --fields \"$vep_fields\""
	eval $cmd >&2
	status=$?
	
	if [ ! $status -eq 0 ]; then
		echo "Error in VEP" >&2
		new_rm $output_fn
		return 1
	fi
	echo -e "`date`: End VEP" >&2
	
	export OMP_NUM_THREADS=1
	[ ! $(which gzip > /dev/null; echo $?) -eq 0 ] && echo "No gzip found" >&2 && return 1
	echo -e "`date`: gzip VEP annotation" >&2
	gzip $output_fn
	new_rm $output_fn
	
	return 0
	
}

srcPL_dnaCallAll=1

## EOF
