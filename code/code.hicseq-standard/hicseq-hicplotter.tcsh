#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-hicplotter.tcsh OUTPUT-DIR PARAM-SCRIPT BRANCH OBJECT(S)
##

if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set objects = ($4)

# test number of input objects
set object = $objects[1]
if ($#objects != 1) then
  scripts-send2err "Error: this operation allows only one input object!"
  exit 1
endif

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir bin_size"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

set inpdir = $branch/$object
set workdir = $outdir/work
mkdir -p $workdir

set chrom_list = `echo $regions | tr ' ' '\n' | cut -d':' -f1 | sort -u | tr '\n' '|' | sed 's/|$//'`
set est_matrices = `cd $inpdir; ls -1 matrix.*.tsv matrix.*.RData | grep -wE "$chrom_list"`            # TODO: inpdir should not have both tsv and RData, check!
foreach est_mat ($est_matrices)
  scripts-send2err "Processing matrix $est_mat..."
  set chr = `echo $est_mat | cut -d'.' -f2`
	
  # extract matrices from RData file
  if (`echo $est_mat | grep -c '\.RData$'` == 1) then
    Rscript ./code/hic-matrix.r matrices -v -o $workdir/tmp $inpdir/$est_mat
  else
    mkdir -p $workdir/tmp
    cat $inpdir/$est_mat >! $workdir/tmp/matrix.k=001.tsv
  endif
  
  # convert matrices and run hicplotter
  foreach mat ($workdir/tmp/matrix.*.tsv)
    set pref = `basename $mat .tsv | sed 's/^matrix\.//'`.$chr
    set inpmat = $pref.matrix.txt
    Rscript $hicplotter_matrix_path/create-hicplotter-matrix.r $workdir/$inpmat $mat 
  end
  rm -rf $workdir/tmp
end

# Run hicplotter
# Go through the regions
foreach region ($regions)
  echo $region
  set chrom = `echo $region | cut -d':' -f1`
  set start = `echo $region | cut -d'-' -f1 | cut -d':' -f2`
  set stop = `echo $region | cut -d'-' -f2`
  set start_bin = `echo $start/$bin_size | bc`
  set stop_bin = `echo $stop/$bin_size | bc`
  set bedgraphs_csv = `echo $bedgraphs | sed 's/ /,/g'`
  set bedgraph_labels_csv = `echo $bedgraph_labels | sed 's/ /,/g'`
  set hic_matrices = `cd $workdir; ls -1 *.matrix.txt | grep -w $chrom | tr '\n' ' '`
  set hic_matrix_no = `echo $hic_matrices | tr ' ' '\n' | wc -l`

  #Get the region labels
  set region_labels = ()
  #Get the loop beds
  set loop_beds = ()
  #Get all bedgraphs
  set all_bedgraphs = ()
  #Get all bedgraph labels
  set all_bedgraph_labels = ()
  #Get histogram type
  set all_histogram_types = ()
  #Get histogram heights
  set all_histogram_heights = ()

  foreach i ( `seq 1 1 $hic_matrix_no` )
    set region_labels = ( $region_labels $region )
    set loop_beds = ( $loop_beds $loop_bed )
    set all_bedgraphs = ( $all_bedgraphs $bedgraphs_csv )
    set all_bedgraph_labels = ( $all_bedgraph_labels $bedgraph_labels_csv )
    set all_histogram_types = ( $all_histogram_types $histogram_type  )
    set all_histogram_heights = ( $all_histogram_heights $histogram_height )
  end

  #Run HiC-Plotter
  set p = `pwd`
  cd $workdir
  python $hicplotter_path/HiCPlotter2.py -v -f $hic_matrices -n $region_labels -chr $chrom -s $start_bin -e $stop_bin -r $bin_size -o $region -hist $all_bedgraphs -hl $all_bedgraph_labels -fhist $all_histogram_types -hm $all_histogram_heights -fh $fileheader -pi $insulation_score #-high $highlight -hf $highlight_bed -peak $loop_beds
  cd $p
end

# Get all the regions in one .pdf
foreach f (`cd $workdir; ls -1 chr*.pdf`)
  mv $workdir/$f $outdir/`echo $f | tr ':' ' ' | tr '-' ' ' | awk -F" " '{print $1":"$2"-"$3}'`.pdf
end

# Set the pdfs in order
set p = `pwd`
cd $outdir
set pdfs = `echo $regions | sed 's/$/ /' | sed 's/ /.pdf /g'` 
#$pdftk_path/pdftk $pdfs cat output plots.pdf
cd $p

# Cleanup
rm -rf $workdir
	
# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


