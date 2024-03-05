#!/bin/sh

[ ! -z $srcPL_genomic ] && [ $srcPL_genomic -eq 1 ] \
	&& return 0

[ -z "$git_dir" ] && git_dir=$(cd $(dirname $BASH_SOURCE)/../..; pwd)

for fn in install linux_perl; do
	. $git_dir/baSHic/scripts/$fn.sh
done

# Installation functions
install_htslib(){
	local version v1 pkg pkg_ver apps_dir status cmd
	local url inst_dir down_dir load_env
	
	install_args $@ -p htslib -d "1.15.1"; status=$?
	[ $status -eq 2 ] && return 0; [ ! $status -eq 0 ] && return 1
	url=https://github.com/samtools/htslib/releases/download
	url=$url/$version/htslib-$version.tar.bz2
	
	# Load environment
	if [ $load_env -eq 1 ]; then
		[ ! -f $inst_dir/lib/pkgconfig/htslib.pc ] \
			&& return 1
		update_env -e PKG_CONFIG_PATH -a "$inst_dir/lib/pkgconfig"
		pkg-config --exists --print-errors htslib >&2 \
			|| return 1
		[ ! -f $inst_dir/bin/htsfile ] \
			&& echo -e "Install $pkg_ver" >&2 \
			&& return 1
		CPPFLAGS="$CPPFLAGS `pkg-config --cflags htslib`"
		LDFLAGS="$LDFLAGS `pkg-config --libs htslib`"
		update_env -e PATH -a "$inst_dir/bin"
		update_env -e LD_LIBRARY_PATH -a "$inst_dir/lib"
		export HTSLIB_DIR=$inst_dir/include/htslib
		return 0
	fi
	
	extract_url -u $url -a $apps_dir -s $pkg_ver
	[ $? -eq 1 ] && return 0
	new_mkdir $inst_dir
	cd $inst_dir
	
	# Set environment
	clear_env
	local CPPFLAGS LDFLAGS
	cmd=$(prep_env_cmd -a $apps_dir -p gcc libtool \
		bzip2 xz zlib curl)
	eval $cmd >&2 || return 1
	
	# Install
	cmd="$down_dir/configure"
	[ ! -z "$CPPFLAGS" ] && cmd="$cmd CPPFLAGS=\"$CPPFLAGS\""
	[ ! -z "$LDFLAGS" ] && cmd="$cmd LDFLAGS=\"$LDFLAGS\""
	cmd="$cmd --prefix=$inst_dir >&2"
	cmd="$cmd && make >&2 && make install >&2"
	eval $cmd
	
	status=$?
	install_wrapup -s $status -i $inst_dir -d $down_dir
	return $status
	
}
install_samtools(){
	local version v1 pkg pkg_ver apps_dir status cmd
	local url inst_dir down_dir load_env
	
	install_args $@ -p samtools -d 1.15.1; status=$?
	[ $status -eq 2 ] && return 0; [ ! $status -eq 0 ] && return 1
	url=https://github.com/samtools/samtools/releases/download
	url=$url/$version/samtools-$version.tar.bz2
	
	# Load environment
	if [ $load_env -eq 1 ]; then
		[ ! -f $inst_dir/bin/samtools ] \
			&& echo -e "Install $pkg_ver" >&2 \
			&& return 1
		update_env -e PATH -a "$inst_dir/bin"
		return 0
	fi
	
	extract_url -u $url -a $apps_dir -s $pkg_ver
	[ $? -eq 1 ] && return 0
	mv $down_dir $inst_dir
	cd $inst_dir
	
	# Set environment
	clear_env
	local CPPFLAGS LDFLAGS
	cmd=$(prep_env_cmd -a $apps_dir -p gcc libtool \
		ncurses xz)
	eval $cmd >&2 || return 1
	
	# Install
	cmd="./configure"
	[ ! -z "$CPPFLAGS" ] && cmd="$cmd CPPFLAGS=\"$CPPFLAGS\""
	[ ! -z "$LDFLAGS" ] && cmd="$cmd LDFLAGS=\"$LDFLAGS\""
	cmd="$cmd --prefix=$inst_dir >&2"
	cmd="$cmd && make >&2 && make install >&2"
	eval $cmd
	
	status=$?
	install_wrapup -s $status -i $inst_dir -d $down_dir
	return $status
	
}
install_bedtools(){
	local version v1 pkg pkg_ver apps_dir status cmd
	local url inst_dir down_dir load_env tmp_dir
	
	install_args $@ -p bedtools -d 2.30.0; status=$?
	[ $status -eq 2 ] && return 0; [ ! $status -eq 0 ] && return 1
	url=https://github.com/arq5x/bedtools2/releases/download
	url=$url/v$version/bedtools-$version.tar.gz
	
	# Load environment
	if [ $load_env -eq 1 ]; then
		[ ! -f $inst_dir/bin/bedtools ] \
			&& echo -e "Install $pkg_ver" >&2 \
			&& return 1
		update_env -e PATH -a "$inst_dir/bin"
		return 0
	fi
	
	extract_url -u $url -a $apps_dir -s $pkg_ver
	[ $? -eq 1 ] && return 0
	tmp_dir=$(ls $apps_dir/downloads | grep $pkg)
	tmp_dir=$apps_dir/downloads/$tmp_dir
	mv $tmp_dir $down_dir; mv $down_dir $inst_dir
	cd $inst_dir
	
	# Set environment
	clear_env
	local CPPFLAGS LDFLAGS
	cmd=$(prep_env_cmd -a $apps_dir -p gcc libtool xz)
	eval $cmd >&2 || return 1
	
	# Install
	make -e VERBOSE=true CPPFLAGS="$CPPFLAGS" \
		LDFLAGS="$LDFLAGS" >&2
	
	status=$?
	install_wrapup -s $status -i $inst_dir -d $down_dir
	return $status
	
}
install_strelka2(){
	local version v1 pkg pkg_ver apps_dir status cmd
	local url inst_dir down_dir load_env tmp_dir
	
	install_args $@ -p strelka -d 2.9.10; status=$?
	[ $status -eq 2 ] && return 0; [ ! $status -eq 0 ] && return 1
	url=https://github.com/Illumina/strelka/releases/download
	url=$url/v$version/strelka-$version.release_src.tar.bz2
	
	# Load environment
	if [ $load_env -eq 1 ]; then
		if [ ! -f $inst_dir/bin/configureStrelkaSomaticWorkflow.py ] \
			|| [ ! -f $inst_dir/bin/configureStrelkaGermlineWorkflow.py ]; then
			echo -e "Install $pkg_ver" >&2 \
				&& return 1
		fi
		update_env -e PATH -a "$inst_dir/bin"
		return 0
	fi
	
	extract_url -u $url -a $apps_dir -s $pkg_ver
	[ $? -eq 1 ] && return 0
	tmp_dir=$(ls $apps_dir/downloads | grep $pkg)
	tmp_dir=$apps_dir/downloads/$tmp_dir
	mv $tmp_dir $down_dir
	new_mkdir $inst_dir
	cd $inst_dir
	
	# Set environment
	clear_env
	local CPPFLAGS LDFLAGS
	cmd=$(prep_env_cmd -a $apps_dir -p gcc libtool \
		zlib cmake boost)
	eval $cmd >&2 || return 1
	
	# Install
	cmd="$down_dir/configure"
	cmd="$cmd --prefix=$inst_dir"
	status=$(which cmake > /dev/null; echo $?)
	[ $status -eq 0 ] && cmd="$cmd --with-cmake=$(which cmake)"
	# cmd="$cmd && make -C $inst_dir >&2"
	cmd="$cmd >&2 && make install >&2"
	eval $cmd
	
	status=$?
	install_wrapup -s $status -i $inst_dir -d $down_dir
	[ ! $status -eq 0 ] && return 1
	
	# Run demos for somatic/germline workflows to double check everything
	echo -e "${cyan}Run Strelka2 test ...${NC}" >&2
	local test_dir=$HOME/strelka_test
	new_rm $test_dir
	new_mkdir $test_dir
	cd $test_dir
	bash $inst_dir/bin/runStrelkaSomaticWorkflowDemo.bash >&2 \
		&& bash $inst_dir/bin/runStrelkaGermlineWorkflowDemo.bash >&2
	status=$?
	if [ $status -eq 0 ]; then
		echo -e "${cyan}Strelka2 test complete${NC}" >&2
		cd
		new_rm $test_dir && return $status
	else
		echo -e "${red}Error with Strelka2 demos!${NC}" >&2 \
			&& return $status
	fi
	
}
install_VEP(){
	# Source: https://m.ensembl.org/info/docs/tools/vep/script/vep_download.html
	local url apps_dir inst_dir perl_dir status
	local release module resp cmd cache_dir
	
	while [ ! -z "$1" ]; do
		case $1 in
			-a | --apps_dir )
				shift
				apps_dir="$1"
				;;
			-c | --cache_dir )
				shift
				cache_dir="$1"
				;;
			-g | --genome )
				shift
				genome="$1"
				;;
			-r | --release )
				shift
				release="$1"
				;;
		esac
		shift
	done
	
	[ -z $apps_dir ] 	&& apps_dir=$HOME/apps
	if [ -z $release ]; then
		make_menu -p "Which release of VEP to install on? (e.g. 105, 106, 111)"
		read release
	fi
	[ -z "$release" ] && echo "Error release missing, exitting" >&2 && return 1
	inst_dir=$apps_dir/vep-$release
	[ -z "$cache_dir" ] && echo "Add -c <cache_dir>" >&2 && return 1
	
	cd $apps_dir
	if [ ! -d $inst_dir ]; then
		git clone https://github.com/Ensembl/ensembl-vep.git >&2
		mv $apps_dir/ensembl-vep $inst_dir
	fi
	
	cd $inst_dir
	git pull >&2
	git checkout release/$release >&2
	
	if [ 1 -eq 2 ]; then
		$perl_dir/bin/perl -I $perl_dir/lib/perl5 -Mlocal::lib=$perl_dir
		eval "$($perl_dir/bin/perl -I $perl_dir/lib/perl5 -Mlocal::lib=$perl_dir)"
		
		module=DBD::mysql
		
		# Remove one module
		$perl_dir/bin/cpanm --uninstall --local-lib=$perl_dir $module
		
		# Install one module
		$perl_dir/bin/cpanm --local-lib=$perl_dir $module
		
		# Check if module successfully installed and location
		$perl_dir/bin/perl -I $perl_dir/lib/perl5 -e "use $module" # check for error
		$perl_dir/bin/perldoc -l $module # location
		
	fi
	
	# Add code to check for mysql_config
	
	# Install Perl modules
	install_perl_modules -a $apps_dir \
		-d expat db bzip2 xz zlib curl htslib \
		-m DBI DBD::mysql Try::Tiny XML::Parser XML::Twig \
		XML::DOM ExtUtils::CBuilder DB_File DB_File::HASHINFO \
		BioPerl Test::Warnings Bio::DB::HTS
		# Archive::Zip
	status=$?
	[ ! $status -eq 0 ] && echo -e "${red}Some perl module failed${NC}" >&2 \
		&& return 1
	
	# Set environment
	clear_env
	local CPPFLAGS LDFLAGS
	cmd=$(prep_env_cmd -a $apps_dir -p gcc libtool perl \
		bzip2 xz zlib curl expat db htslib)
	eval $cmd >&2 || return 1
	
	# VEP install and download cached database files (# 460) and plugins (gnomad)
	[ ! -d "$cache_dir" ] && mkdir "$cache_dir"
	cd $inst_dir
	cmd="perl INSTALL.pl --NO_HTSLIB --CACHEDIR $cache_dir"
	eval $cmd >&2
	[ ! $? -eq 0 ] && echo -e "Error in VEP installation" >&2 && return 1
	
	return 0
	
}

srcPL_genomic=1

###

