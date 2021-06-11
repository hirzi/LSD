#!/bin/bash

##### This is a bash-based msToGLF and ANGSD wrapper that simulates genotype likelihoods and calculates summary statistics from msms output. It has been written to be compatible with ABCtoolbox.
##### Hirzi Luqman, 11.03.2019
##### Example usage: ./lsd_low.sh -f neutral_demo_simpleModel_11params_RECON2.out -p 20,20,20,20,20,20 -l 5000 -d 2 -e 0.01 -r ${REF_index} -w ${working_dir} -o results.concatenated
##### Recall, -p takes diploid sample size (not haploid!)

##### For reference, see:
#http://popgen.dk/angsd/index.php/MsToGlf
#http://popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests
#https://github.com/mfumagalli/Tjarno/blob/master/Files/selection_2.md


##### Load required modules
module load gcc/4.9.2 gdc angsd/0.925


##### Parse input arguments
while getopts f:p:l:d:e:r:w:o: option
do
case "${option}"
in
f) ms_file=${OPTARG};;
p) pop_info=${OPTARG};;
l) sequence_length=${OPTARG};;
d) depth=$OPTARG;;
e) error_rate=$OPTARG;;
r) REF_index=$OPTARG;;
w) working_dir=$OPTARG;;
o) output_name=$OPTARG;;
esac
done

if [ ! -d ${working_dir} ]; then
	mkdir ${working_dir}
fi


##### Process input variables 
# Name prefix
prefix=$(basename "$ms_file")
# Split input population string into population array (https://stackoverflow.com/questions/10586153/split-string-into-an-array-in-bash)
IFS=', ' read -r -a pop_array <<< "$pop_info"
# Count number of total individuals, aka sum of array (https://stackoverflow.com/questions/13635293/unix-shell-script-adding-the-elements-of-an-array-together)
no_inds_total=0
for i in ${pop_array[@]}; do
  let no_inds_total+=$i
done
# Number of populations, aka length of array
no_pops=${#pop_array[@]}

# Print input (verbose)
echo "Input file:" $ms_file
echo "Reference index is given at:" $REF_index
echo "The working directory is given as:" $working_dir
echo "The final output file will be named:" $output_name
echo "Population array:" $pop_info
echo "Simulated sequence length:" $sequence_length"bp"
echo "Simulated depth:" $depth"X"
echo "Simulated error rate:" $error_rate
echo "Total number of individuals:" $no_inds_total
echo "Total number of populations:" $no_pops


##### Generate a list of populations, and of all population pairs.
if [ ! -f ${working_dir}/pop_list ]; then
	# First make a list of population names
	echo pop1 > ${working_dir}/pop_list
	for pop in $(seq 2 ${no_pops}); do
		echo pop${pop} >> ${working_dir}/pop_list
	done

	# Then make pairwise list from the list above:
	set -- $(cat ${working_dir}/pop_list | cut -d$'\t' -f1)
	for a; do
	    shift
	    for b; do
	        printf "%s\t%s\n" "$a" "$b"
	#		printf "$a - $b\n"
	    done
	done > ${working_dir}/pop_name_pairs
fi
# Number of population pairs
num_pairs=$(cat ${working_dir}/pop_name_pairs | wc -l)


##### Generate null output template, in the case of a null (zero segsites) or invalid (segsites > sequence length) msms output
if [ ! -f ${working_dir}/null.headers ]; then
	for pop in $(seq 1 ${no_pops}); do
		# To add suffix after each word, see: https://stackoverflow.com/questions/28984295/add-comma-after-each-word
		echo -e "tW"'\t'"tP"'\t'"Tajima" | sed "s/\>/.pop${pop}/g" > ${working_dir}/pop${pop}.null_headers.thetas
		echo -e "singletons"'\t'"doubletons" | sed "s/\>/.pop${pop}/g" > ${working_dir}/pop${pop}.null_headers.sfs
	done
	paste ${working_dir}/pop*.null_headers.thetas > ${working_dir}/null_headers.thetas
	paste ${working_dir}/pop*.null_headers.sfs > ${working_dir}/null_headers.sfs

	for pop in $(seq 1 ${num_pairs}); do
		pop_pair=`sed -n ${pop}p < ${working_dir}/pop_name_pairs`
		# Remember, set allows you to define the elements of your list as variables, according to their order
		set -- $pop_pair
		# Print headers
		echo "${1}.${2}.fst" > ${working_dir}/${1}.${2}.null_headers.fst
	done
	paste ${working_dir}/pop*.pop*.null_headers.fst > ${working_dir}/null_headers.fst
	# Concatenate all headers
	paste ${working_dir}/null_headers.thetas ${working_dir}/null_headers.sfs ${working_dir}/null_headers.fst > ${working_dir}/null.headers
	# Remove temporary files
	rm ${working_dir}/pop*.null_headers.* ${working_dir}/null_headers.*
fi
# Number of headers
num_headers=$(cat ${working_dir}/null.headers | awk '{print NF}')


##### MAIN FUNCTION #####
## Condition flow control and output according to whether msms output contains non-zero segsites and comprises number of segsites < sequence length. If not, write out null output files.
# If the msms output has no segsites, then it will have less than 5 lines.
nosegsites_linecount_threshold=5
len_msout=$(cat ${working_dir}/${ms_file} | wc -l)
num_segsites=$(cat ${working_dir}/${ms_file} | sed -n 5p | awk '{print $2}')

# If the msms output has no segsites, then it will have less than 5 lines. We generate an all zeroes output in this case.
if [ "$len_msout" -lt "$nosegsites_linecount_threshold" ]; then
	echo "msms simulation produced zero segregating sites - generating null output file (all 0s)"
	cat ${working_dir}/null.headers > ${working_dir}/${output_name}
	printf '0\t%.0s' $(seq ${num_headers}) >> ${working_dir}/${output_name}

# If the msms output has segsites > sequence length, the result is invalid (and can't be read by msToGlf). Here, the output is given all 999999 (i.e. we designate this as the NA value, which we can remove/filter out later)
elif [ "$len_msout" -ge "$nosegsites_linecount_threshold" ] && [ "$num_segsites" -ge "$sequence_length" ]; then
	echo "msms simulation produced invalid number of segregating sites (more segregating sites than designated sequence length) - generating null output file (all NAs/999999s)"
	cat ${working_dir}/null.headers > ${working_dir}/${output_name}
	printf '999999\t%.0s' $(seq ${num_headers}) >> ${working_dir}/${output_name}

# If the msms output has non-zero segsites, and has segsites < sequence_length, we proceed as follows:
elif [ "$len_msout" -ge "$nosegsites_linecount_threshold" ] && [ "$num_segsites" -lt "$sequence_length" ]; then
	echo "msms simulation produced non-zero and valid number of segregating sites - proceeding with analysis..."
	
	##### Simulate genotype likelihoods (msToGLF)
	# Simulate genotype likelihoods from msms output via msToGlf defining error rate, sequencing depth (you can also supply a depthFile filename if you want to force a different mean depth between individuals), number of diploid individuals, and sequence length (simulating invariant sites). Outputting a single replicate from same scenario (-singleOut 1)
	echo "Simulating genotype likelihoods with error rate of" $error_rate "and depth of" $depth"X; also simulating invariant sites assuming sequence length of" $sequence_length"bp"
	msToGlf -in ${working_dir}/${ms_file} -out ${working_dir}/${prefix}.gl -err ${error_rate} -depth ${depth} -nind ${no_inds_total} -singleOut 1 -regLen ${sequence_length} &> /dev/null

	##### Calculate single population summary statistics (ANGSD)
	end=0
	for pop in $(seq 1 ${no_pops}); do
		echo "Processing population" $pop
		# We define start and end to split the simulated glf according number of individuals per population. 
		let end+="${pop_array[${pop}-1]}"
		start=$(( ${end} - ${pop_array[${pop}-1]} + 1 ))
		no_inds_per_pop="${pop_array[${pop}-1]}"
		echo "Population" $pop "- starting individual:" $start "; ending individual:" $end "; population size:" $no_inds_per_pop
		# Extract populations from GLF file (by number of individuals per population)
		splitgl ${working_dir}/${prefix}.gl.glf.gz ${no_inds_total} ${start} ${end} > ${working_dir}/pop${pop}.glf.gz &> /dev/null
		# Calculate the site allele frequency likelihoods (doSaf)
		echo "Calculating SAF for population" $pop
		angsd -glf ${working_dir}/pop${pop}.glf.gz -nInd ${no_inds_per_pop} -doSaf 1 -isSim 1 -P 1 -fai ${REF_index} -out ${working_dir}/pop${pop} &> /dev/null
		#angsd -glf ${working_dir}/pop${pop}.glf.gz -nInd ${no_inds_per_pop} -doSaf 1 -doMajorMinor 1 -doMaf 1 -isSim 1 -P 1 -fai ${REF_index} -out ${working_dir}/pop${pop}
		# Calculate the SFS (for use as prior in calculation of thetas)
		echo "Calculating SFS for population" $pop
		realSFS ${working_dir}/pop${pop}.saf.idx -P 1 2> dev/null > ${working_dir}/pop${pop}.sfs 
		# Calculate the site allele frequency likelihoods (doSaf)
		echo "Calculating thetas for population" $pop	
		# angsd -glf ${working_dir}/pop${pop}.glf.gz -nInd ${no_inds_per_pop} -doSaf 1 -doThetas 1 -isSim 1 -P 1 -pest ${working_dir}/pop${pop}.sfs -fai ${REF_index} -out ${working_dir}/pop${pop}
		realSFS saf2theta ${working_dir}/pop${pop}.saf.idx -P 1 -sfs ${working_dir}/pop${pop}.sfs -outname ${working_dir}/pop${pop} &> /dev/null
		#angsd -glf ${working_dir}/pop${pop}.glf.gz -nInd ${no_inds_per_pop} -doSaf 1 -doMajorMinor 1 -doMaf 1 -doThetas 1 -isSim 1 -P 1 -pest ${working_dir}/pop${pop}.sfs -fai ${REF_index} -out ${working_dir}/pop${pop}
		thetaStat do_stat ${working_dir}/pop${pop}.thetas.idx &> /dev/null
	done
	echo "All populations' thetas calculated!"


	# ##### Calculate pairwise population summary statistics (ANGSD)
	# echo "Preparing calculation of pairwise statistics..."
	# for pop in $(seq 1 ${num_pairs}); do
	# 	pop_pair=`sed -n ${pop}p < ${working_dir}/pop_name_pairs`
	# 	echo "Processing population pair" $pop_pair
	# 	# Remember, set allows you to define the elements of your list as variables, according to their order
	# 	set -- $pop_pair
	# 	# Calculate the 2DSFS prior
	# 	echo "Calculating 2D SFS for population pair" $pop_pair	
	# 	realSFS ${working_dir}/${1}.saf.idx ${working_dir}/${2}.saf.idx -P 1 > ${working_dir}/${1}.${2}.ml
	# 	# Calculate the FST
	# 	echo "Calculating FST for population pair" $pop_pair
	# 	realSFS fst index ${working_dir}/${1}.saf.idx ${working_dir}/${2}.saf.idx -sfs ${working_dir}/${1}.${2}.ml -fstout ${working_dir}/${1}.${2}.stats -whichFst 1
	# 	# Get the global estimate (here we output only the weighted estimate)
	# 	echo "${1}.${2}.fst" > ${working_dir}/${1}.${2}.globalFST
	# 	realSFS fst stats ${working_dir}/${1}.${2}.stats.fst.idx 2> /dev/null | cut -f 2 >> ${working_dir}/${1}.${2}.globalFST
	# 	# Get sliding window estimates
	# 	#realSFS fst stats2 ${working_dir}/${1}.${2}.stats.fst.idx -win 500 -step 500 > ${OUT}/FST_results/slidingwindow_${1}.${2}
	# done
	# echo "All populations' FSTs calculated!"

	echo "Preparing calculation of pairwise statistics..."
		for pop in $(seq 1 ${num_pairs}); do
			pop_pair=`sed -n ${pop}p < ${working_dir}/pop_name_pairs`
			echo "Processing population pair" $pop_pair
			# Remember, set allows you to define the elements of your list as variables, according to their order
			set -- $pop_pair
			echo "${1}.${2}.fst" > ${working_dir}/${1}.${2}.globalFST
			# Calculate the 2DSFS prior
			echo "Calculating 2D SFS for population pair" $pop_pair	
			realSFS ${working_dir}/${1}.saf.idx ${working_dir}/${2}.saf.idx -P 1 2> /dev/null > ${working_dir}/${1}.${2}.ml 
			if [[ ! ${no_pops} = 3 ]]; then
				echo "there are not 3 but" ${no_pops} "populations, won't do PBS"
				# Calculate the FST
				 echo "Calculating FST for population pair" $pop_pair
				realSFS fst index ${working_dir}/${1}.saf.idx ${working_dir}/${2}.saf.idx -sfs ${working_dir}/${1}.${2}.ml -fstout ${working_dir}/${1}.${2}.stats -whichFst 1 &> /dev/null
				# Get the global estimate (here we output only the weighted estimate)
				realSFS fst stats ${working_dir}/${1}.${2}.stats.fst.idx 2> /dev/null | cut -f 2 >> ${working_dir}/${1}.${2}.globalFST
				 echo "All populations' FSTs calculated!"
				# Get sliding window estimates
				#realSFS fst stats2 ${working_dir}/${1}.${2}.stats.fst.idx -win 500 -step 500 > ${OUT}/FST_results/slidingwindow_${1}.${2}
			fi
		done


	###### Calculate 3-population summary statistics, e.g. PBS (ANGSD) - ONLY IF THERE'S A RELEVANT OUTGROUP & A SUITABLE SAMPLE DESIGN - hashed out for now (needs a little revision)
	## Revision note 1: Need to add collation script for PBS below.
	## Revision note 2: If implemented, make sure to hash out calculation of pairwise summary statistics above, as the PBS part already calculates FSTs. Streamline and code accordingly. 
	# if [[ ${no_pops} > 2 ]]; then
	# 	# First, generate a list of population trios. Note this generates a list of all population trios (combinations). In reality, you'll want to have a single outgroup population for which to perform PBS, hence you'll rather append the outgroup population to each population pair in pop_name_pairs.
	# 	# Revision note 3: Don't forget to adjust python code below to accommodate arbitrary number of populations, and to specify custom working_dir.
	# 	# python - <<-EOF
	# 	# 	import itertools
	# 	# 	trios=list(itertools.combinations(["pop1","pop2","pop3","pop4","pop5","pop6"],3))
	# 	# 	with open('pop_name_trios', 'w') as fp:
	# 	# 	    fp.write('\n'.join('{}\t{}\t{}'.format(x[0],x[1], x[2]) for x in trios))
	# 	# EOF
	# 	# Append final newline character
	# 	sed -i -e '$a\' ${working_dir}/pop_name_trios
	# 	# Number of population trios
	# 	num_trios=$(cat ${working_dir}/pop_name_trios | wc -l)
	# 	# Then calculate FST and PBS (recall ANGSD calculates PBS from FSTs, hence it calculates these statistics together if 3 populations are supplied)
	# 	for pop in $(seq 1 ${num_trios}); do
	# 		pop_trio=`sed -n ${pop}p < ${working_dir}/pop_name_trios`
	# 		echo "Processing population trio" $pop_trio
	# 		# Remember, set allows you to define the elements of your list as variables, according to their order
	# 		set -- $pop_trio
	# 		##calculate pbs and fst
	# 		echo "Calculating PBS for population trio" $pop_trio
	# 		realSFS fst index ${working_dir}/${1}.saf.idx ${working_dir}/${2}.saf.idx ${working_dir}/${3}.saf.idx -sfs ${working_dir}/${1}.${2}.ml -sfs ${working_dir}/${1}.${3}.ml -sfs ${working_dir}/${2}.${3}.ml -fstout ${working_dir}/${1}.${2}.${3}.stats -whichFst 1

	# 		#get the global estimate (here we're interested in outputting the PBS results)
	# 		realSFS fst stats ${working_dir}/${1}.${2}.${3}.stats.fst.idx >> ${working_dir}/${1}.${2}.${3}.globalPBS
	# 		echo
	# 	done
	# fi

	if [[ ${no_pops} = 3 ]]; then
	echo "Running in 3-population mode"
		# First, generate a list of population trios. Note this generates a list of all population trios (combinations). In reality, you'll want to have a single outgroup population for which to perform PBS, hence you'll rather append the outgroup population to each population pair in pop_name_pairs.
		# Revision note 3: Don't forget to adjust python code below to accommodate arbitrary number of populations, and to specify custom working_dir.
		# python - <<-EOF
			# import itertools
			# trios=list(itertools.combinations(["pop1","pop2","pop3","pop4","pop5","pop6"],3))
			# with open('pop_name_trios', 'w') as fp: fp.write('\n'.join('{}\t{}\t{}'.format(x[0],x[1], x[2]) for x in trios))
		# EOF
		# Append final newline character
		printf "pop1\tpop2\tpop3" > ${working_dir}/pop_name_trios
		# sed -i -e '$a\' ${working_dir}/pop_name_trios
		# Number of population trios
		# num_trios=$(cat ${working_dir}/pop_name_trios | wc -l)
		# Then calculate FST and PBS (recall ANGSD calculates PBS from FSTs, hence it calculates these statistics together if 3 populations are supplied)
		# for pop in $(seq 1 ${num_trios}); do
		pop=1
			pop_trio=`sed -n ${pop}p < ${working_dir}/pop_name_trios`
			echo "Processing population trio" $pop_trio
			# Remember, set allows you to define the elements of your list as variables, according to their order
			set -- $pop_trio
			## calculate 2DSFS priors
			# echo "calculating 2DSFS for every pair"
			# realSFS ${working_dir}/${1}.saf.idx ${working_dir}/${2}.saf.idx -P 1  2> /dev/null > ${working_dir}/${1}.${2}.ml
			# realSFS ${working_dir}/${1}.saf.idx ${working_dir}/${3}.saf.idx -P 1  2> /dev/null > ${working_dir}/${1}.${3}.ml
			# realSFS ${working_dir}/${2}.saf.idx ${working_dir}/${3}.saf.idx -P 1  2> /dev/null > ${working_dir}/${2}.${3}.ml
			##calculate pbs and fst
			echo "Calculating pairwise FSTs and PBS for population trio" 
			realSFS fst index ${working_dir}/${1}.saf.idx ${working_dir}/${2}.saf.idx ${working_dir}/${3}.saf.idx -sfs ${working_dir}/${1}.${2}.ml -sfs ${working_dir}/${1}.${3}.ml -sfs ${working_dir}/${2}.${3}.ml -fstout ${working_dir}/${1}.${2}.${3}.stats -whichFst 1 &> /dev/null
			#get the global estimate
			realSFS fst stats ${working_dir}/${1}.${2}.${3}.stats.fst.idx >> ${working_dir}/${1}.${2}.${3}.globalFST_PBS
			echo
			head -n1 ${working_dir}/${1}.${2}.${3}.globalFST_PBS | cut -f 2 >> ${working_dir}/${1}.${2}.globalFST
			head -n2 ${working_dir}/${1}.${2}.${3}.globalFST_PBS | tail -n1 | cut -f 2 >> ${working_dir}/${1}.${3}.globalFST
			head -n3 ${working_dir}/${1}.${2}.${3}.globalFST_PBS | tail -n1 | cut -f 2 >> ${working_dir}/${2}.${3}.globalFST
		# done
	fi

	##### Collate results
	echo "Collating results and generating final output..."

	# Collate theta results
	for pop in $(seq 1 ${no_pops}); do	
		pop_suffix=".pop${pop}"
		# For adding suffix, see: https://unix.stackexchange.com/questions/265335/adding-a-number-as-a-suffix-to-multiple-columns. For expanding variable in awk, use the -v option; see: https://unix.stackexchange.com/questions/340369/expanding-variables-in-awk. Tab separate fields via -v OFS='\t' (see: https://askubuntu.com/questions/231995/how-to-separate-fields-with-space-or-tab-in-awk)
		cat ${working_dir}/pop${pop}.thetas.idx.pestPG | head -n 1 | awk  -v OFS='\t' '{ print $4, $5, $9 }' | awk -v OFS='\t' -v pop_suffix="$pop_suffix" '{$1 = $1 pop_suffix; $2 = $2 pop_suffix; $3 = $3 pop_suffix; print }' > ${working_dir}/pop${pop}.thetas.idx.headers
		# we need to get mean values across loci, whereas lsd_low only uses one (simulated) locus so tail is point estimates		
		# cat ${working_dir}/pop${pop}.thetas.idx.pestPG | tail -n 1 | awk  -v OFS='\t' '{ print $4, $5, $9 }' >> ${working_dir}/pop${pop}.thetas.idx.headers
		mean_tW=$(tail -n+1 ${working_dir}/pop${pop}.thetas.idx.pestPG | awk -v OFS='\t' '{ sum += $4 } END { if (NR > 0) print sum / NR }' )
		mean_tP=$(tail -n+1 ${working_dir}/pop${pop}.thetas.idx.pestPG | awk -v OFS='\t' '{ sum += $5 } END { if (NR > 0) print sum / NR }' )
		mean_TajimaD=$(tail -n+1 ${working_dir}/pop${pop}.thetas.idx.pestPG | awk -v OFS='\t' '{ sum += $9 } END { if (NR > 0) print sum / NR }' )
		printf "$mean_tW\t$mean_tP\t$mean_TajimaD" >> ${working_dir}/pop${pop}.thetas.idx.headers
	done
	paste ${working_dir}/pop*.thetas.idx.headers > ${working_dir}/thetas.results.concatenated

	# Collate SFS results (here retaining just singleton and doubleton categories)
	for pop in $(seq 1 ${no_pops}); do
		echo -e "singletons.pop"${pop}'\t'"doubletons.pop"${pop} > ${working_dir}/pop${pop}.SFS.truncated.results
		cat ${working_dir}/pop${pop}.sfs | awk -v OFS='\t' '{ print $2, $3 }' >> ${working_dir}/pop${pop}.SFS.truncated.results
	done
	paste ${working_dir}/pop*.SFS.truncated.results > ${working_dir}/SFS.truncated.results.concatenated

	# Collate FST results
	paste ${working_dir}/pop*.pop*.globalFST > ${working_dir}/fst.results.concatenated

	# Collate PBS results
	if [ -f ${working_dir}/${1}.${2}.${3}.globalFST_PBS ];then
	tail -n+4 ${working_dir}/${1}.${2}.${3}.globalFST_PBS | \
	awk OFS='\t' '{for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }} NF>p { p = NF }
	END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
		}' | sed 's/ /\t/g' > ${working_dir}/PBS.results.concatenated
		# Concatenate all results (make sure all fields are tab-separated)
	paste ${working_dir}/thetas.results.concatenated ${working_dir}/SFS.truncated.results.concatenated ${working_dir}/fst.results.concatenated ${working_dir}/PBS.results.concatenated > ${working_dir}/${output_name}
	else
		# Concatenate all results (make sure all fields are tab-separated)
	paste ${working_dir}/thetas.results.concatenated ${working_dir}/SFS.truncated.results.concatenated ${working_dir}/fst.results.concatenated > ${working_dir}/${output_name}
	fi


	# Check that number of fields match between headers and results
	#cat ${working_dir}/${output_name} | awk '{print NF}'
	echo "GLF simulation and calculation of summary statistics has finished!"

	##### Remove temporary and intermediate files
	#rm

# For all other cases (there shouldn't be any, but just in case), we output an all 999999 (NA) file.
else
	echo "msms simulation produced invalid number of segregating sites - generating null output file (all NAs/999999s)"
	cat ${working_dir}/null.headers > ${working_dir}/${output_name}
	printf '999999\t%.0s' $(seq ${num_headers}) >> ${working_dir}/${output_name}
fi

##### END OF MAIN FUNCTION #####

