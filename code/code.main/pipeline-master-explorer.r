#!/usr/bin/Rscript
#$ -S /usr/bin/Rscript

########################################################################################################################################################
usage = 'pipeline-master-explorer.r [OPTIONS] SCRIPT OUTDIR-PREFIX PARAM-SCRIPTS INPUT-BRANCHES SPLIT-VARIABLE OUTPUT-OBJECT-VARIABLE TUPLES' 
########################################################################################################################################################



## 
## check_status
##
check_status = function(v)
{
  missing = min(file.exists(v["inp-branch"]))==0
  if (missing==TRUE) return("missing-inputs")
  inp_files = paste(v["inp-branch"],unlist(strsplit(v["inp-object"],split="[, ]")),sep='/')
  if (v["inp-object"]=="") inp_files = Sys.glob(paste(inp_files,"*",sep="/"))
  inp_files = c(inp_files, v["params"])
  missing = min(file.exists(inp_files))==0
  if (missing==TRUE) return("missing-inputs")
  exists = file.exists(v["out-dir"])
  if (exists==FALSE) return("new")
  uptodate = file.info(paste(v["out-dir"],"job.sh",sep='/'))$mtime > max(file.info(inp_files)$mtime)
  if (uptodate==FALSE) return("out-of-date")
  return("up-to-date")
}


## 
## create_out_branch
##
create_out_branch = function(v,ignore_group)
{
  params_name = sub('.[^.]+$','',sub('^params/params.','',v[,"params"]))                          # remove: (a) "params/params." prefix (b) any suffix
  branch_name = sub('^.*/results[/]*','',v[,"inp-branch"])                                        # remove from beginning to "results/"            
  prefix = v[,"out-prefix"]                                                                       # simple prefix
  if (ignore_group==FALSE) prefix = paste(prefix,paste("by_",v[,"out-obj-var"],sep=''),sep='.')   # full prefix (output object variable is added)
  out = paste(paste(prefix,params_name,sep='.'),branch_name,v[,"split-value"],sep='/')
  out = sub('[/]+$','',out)                                                                       # remove trailing slashes
  out = sub("[/]{2,}","/",out)                                                                    # suppress repetitive slashes
  return(out)
}


## 
## create_out_branch_short
##
create_out_branch_short = function(v)
{
  branch_name = sub('^[^/]+/results[/]*','',v$"inp-branch")
  prefix = paste("by_",v$"out-obj-var",sep='')
  out = paste(prefix,branch_name,v$"split-value",sep='/')
  out = sub('[/]+$','',out)                                                                     # remove trailing slashes
  out = gsub("[/]{2,}","/",out)                                                                 # suppress all (gsub) repetitive slashes
  return(out)
}


## 
## create_out_dir
##
create_out_dir <- function(v)
{
  return(paste(v[,"out-branch"],v[,"out-object"],sep='/'))
}


## 
## generate_command
##
generate_command <- function(v)
{
  out_dir = paste(v[,"out-branch"],v[,"out-object"],sep='/')
  cmd = paste(v[,"script"],out_dir,v[,"params"],v[,"inp-branch"],v[,"inp-object"],sep=' ')
  return(cmd)
}



## 
## combn_obj
##
combn_obj = function(mat,tuples) {
  mat1 = aggregate(mat,by=list(mat[,"out-object"]),FUN=function(x) paste(unique(x),collapse=','))[-1]                                 # aggregate entries by output-object
  if (tuples>0) {
    inp_obj_grid = as.matrix(expand.grid(rep(list(mat1[,"inp-object"]),tuples)))                                                        # (input-objects)^N
    out_obj_grid = as.matrix(expand.grid(rep(list(mat1[,"out-object"]),tuples)))                                                        # (output-objects)^N
    obj_mapping = cbind("out-object"=apply(out_obj_grid,1,paste,collapse='.'),"inp-object"=apply(inp_obj_grid,1,paste,collapse=' '))    # input-output object mapping
    out = obj_mapping
    for (k in setdiff(colnames(mat1),c("inp-object","out-object","out-branch-short"))) { out = cbind(out,mat1[1,k]); colnames(out)[ncol(out)] = k }
  } else {
    out = mat1[,colnames(mat1)!="out-branch-short"]
  }
  return(out)
}



##
## create_output_objects
##
create_output_objects = function(inp_db, sample_sheet, out_obj_var, split_var)
{
  # split obj-db by inp-obj-var
  if (opt$verbose) write(paste("Splitting objects by '",split_var,"' and grouping by '",out_obj_var,"'...",sep=''),stderr())
  obj_db = as.data.frame(cbind(inp_db,"out-obj-var"=out_obj_var,"split-var"=split_var))
  if (split_var=="") obj_db = cbind(obj_db,"split-value"="")
  L = split(obj_db,obj_db$"inp-obj-var")                                                                              # separate matrices of different inp-obj-var
  f1 = function(a) { merge(L[[a]],cbind("inp-object"=sample_sheet[,a],"out-object"=sample_sheet[,out_obj_var])) }
  f2 = function(a) { merge(L[[a]],cbind("inp-object"=sample_sheet[,a],"split-value"=sample_sheet[,split_var])) }
  for (a in names(L)) {
    L[[a]] = f1(a)                                                                                                    # merge obj matrix with inp-obj-var/out-obj-var columns from sample sheet
    if (split_var!="") L[[a]] = f2(a)                                                                                 # merge obj matrix with inp-obj-var/split-var columns from sample sheet
  }
  obj_db = do.call(rbind,L)

  # generate input-output object maps
  if (opt$verbose) write("Generating input-output object maps...",stderr())
  obj_db = cbind(obj_db,"out-branch-short"=create_out_branch_short(obj_db))                                           # add output branch (short version)
  L = split(obj_db,obj_db$"out-branch-short")                                                                         # separate entries that have different short output branch

  # run combn on list elements and collect results into a single table
  obj_db = do.call(rbind,lapply(L,combn_obj,tuples=ifelse(out_obj_var=="*",0,tuples)))
  
  return(obj_db)
}



## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## MAIN
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# install packages
for (p in c("optparse")) 
  if (!suppressPackageStartupMessages(require(p,character.only=TRUE,quietly=TRUE,warn.conflicts=FALSE))) {
    install.packages(p,repos="http://cran.rstudio.com/") 
    library(p,character.only=TRUE,quietly=TRUE,verbose=FALSE)
  }

# process command-line arguments
option_list <- list(
  make_option(c("-v","--verbose"), action="store_true",default=FALSE, help="Print more messages."),
  make_option(c("-S","--sample-sheet"), default="inputs/sample-sheet.tsv", help="Sample sheet file name (required) [default \"%default\"]."),
  make_option(c("-F","--filter-branch"), default="", help="Regular expression for filtering input branches [default \"%default\"]."),
  make_option(c("--exclude-branch"), default="", help="Regular expression for excluding input branches [default \"%default\"]."),
  make_option(c("--exclude-obj"), default="", help="Regular expression for excluding input objects [default \"%default\"].")
)
  
# get command line options (if help option encountered print help and exit)
arguments = parse_args(args=commandArgs(trailingOnly=T), OptionParser(usage=usage,option_list=option_list), positional_arguments=c(0,Inf))
opt = arguments$options
inputs = arguments$args
if (length(inputs)!=7) { write(paste("USAGE: ",usage,sep=''),stderr()); quit(save='no') }
script = inputs[1]
out_prefix = inputs[2]
params = Sys.glob(unlist(strsplit(inputs[3],split=' ')))
inp_dirs = Sys.glob(unlist(strsplit(inputs[4],split=' ')))
split_var = inputs[5]
out_obj_vars = inputs[6]
tuples = as.integer(inputs[7])
out_dir = paste(sub("/[^/]+$","",out_prefix),".db",sep='/')

# remove .db directory if it exists, then start from scratch
if (file.exists(out_dir)) { write("Warning: .db directory exists, removing!",stderr()); unlink(out_dir,recursive=TRUE) }
dir.create(out_dir)

# error checking
if (file.exists(opt$"sample-sheet")==FALSE) { write("Error: sample sheet file does not exist!",stderr()); quit(save='no') }
if (length(params)==0) { write("Error: parameter file list is empty!",stderr()); quit(save='no') }


# read and process sample sheet
if (opt$verbose) write("Reading sample sheet...",stderr())
sample_sheet = read.table(opt$"sample-sheet",header=T,check.names=FALSE,as.is=TRUE)
sample_sheet = as.matrix(cbind('*'=rep("all-samples",nrow(sample_sheet)),sample_sheet))                                                 # add a universal identifier to sample sheet

# determine input branches
if (opt$verbose) write("Mapping input/output objects...",stderr())
inp_db = c()
for (inp_dir in inp_dirs) { 
  inp_db_file = paste(inp_dir,"results/.db/db.RData",sep="/")
  if (file.exists(inp_db_file)==FALSE) {
    obj = expand.grid("inp-branch"=paste(inp_dir,'results',sep='/'), "inp-obj-var"="sample", "inp-object"=sample_sheet[,"sample"])      # identify input objects
  } else {
    e = new.env()
    load(inp_db_file,e)
    if (length(e$obj_db)==0) { obj = c()
    } else {
      if ((e$tuples!=1)&&(out_obj_vars!="*")&&(out_obj_vars!=".")) { write("Error: output object variable can only be set to '*' or '.' for combinatorial input objects!", stderr()); quit(save='no') }
      obj = e$obj_db[,c("out-branch","out-obj-var","out-object")]
      colnames(obj) = c("inp-branch","inp-obj-var","inp-object")                                                                        # output objects of previous step now become input objects
      obj[,"inp-branch"] = paste(inp_dir,obj[,"inp-branch"],sep='/')                                                                    # add input directory prefix to input branch
    }
  }
  inp_db = rbind(inp_db,obj)
}

# filter input branches
if (length(inp_db)==0) { obj_db = c()
} else {
  if (opt$verbose) write("Filtering input branches/objects...",stderr())
  selected = grep(opt$"filter-branch",inp_db$"inp-branch",ignore.case=TRUE)
  if (opt$"exclude-branch"!="") selected = intersect(selected, grep(opt$"exclude-branch",inp_db$"inp-branch",ignore.case=TRUE,invert=TRUE))
  if (opt$"exclude-obj"!="") selected = intersect(selected, grep(opt$"exclude-obj",inp_db$"inp-object",ignore.case=TRUE,invert=TRUE))
  if (length(selected)>0) { obj_db = inp_db[selected,,drop=FALSE] } else { obj_db = c() }
}

# check if empty
if (length(obj_db)==0) {
  obj_db = c()
  write("Warning: no input branches found (after filtering), creating empty tables!",stderr())
  save.image(file=paste(out_dir,"db.RData",sep="/"))
  write.table(c(),col.names=FALSE,row.names=FALSE,sep='\t',quote=FALSE,file=paste(out_dir,"run",sep="/"))
  write.table(c(),col.names=FALSE,row.names=FALSE,sep='\t',quote=FALSE,file=paste(out_dir,"db.tsv",sep="/"))
  if (opt$'verbose'==TRUE) write('Done.',stderr())
  quit(save='no')
}

# print output object matrix
#if (opt$verbose) { write("INPUT OBJECT TABLE: ",stderr()); write.table(obj_db,quote=FALSE,col.names=TRUE,row.names=TRUE,sep='\t',stderr()) }

# check if out-object-variable is "*" (i.e. all input objects are grouped)
if ((out_obj_vars=="*")&&(split_var=="")) {
  obj_db[,"inp-obj-var"] = "."
  obj_db[,"inp-object"] = ""
  obj_db = unique(obj_db)
  obj_db = cbind(obj_db,"out-obj-var"="*","split-var"="","split-value"="","out-object"="all-samples")

# check if out-object-variable is "." (i.e. output objects are the same as input objects)
} else if (out_obj_vars==".") {
  if (split_var!="") write("Warning: split-var is ignored!",stderr())                           # TODO: create separate case for split-var!=""
  obj_db = cbind(obj_db,"out-obj-var"=obj_db[,"inp-obj-var"],"split-var"="","split-value"="","out-object"=obj_db[,"inp-object"])

} else {
  L = unlist(strsplit(out_obj_vars,split=' '))
  obj_db = do.call(rbind, lapply(L, function(out_obj_var) create_output_objects(inp_db=obj_db,sample_sheet=sample_sheet,out_obj_var=out_obj_var,split_var=split_var)))
}

# add parameters
if (opt$verbose) write("Adding parameter information...",stderr())
obj_db = do.call(rbind,lapply(params,function(p) cbind(obj_db,"params"=p)))

# add script and prefix information
if (opt$verbose) write("Adding script and prefix...",stderr())
obj_db = cbind(obj_db,"script"=script,"out-prefix"=out_prefix)

# generate output branches and directories
if (opt$verbose) write("Generating output branches and directories...",stderr())
ignore_group = out_obj_vars=="*"
obj_db = cbind(obj_db,"out-branch"=create_out_branch(obj_db,ignore_group=ignore_group))                                               # add output branch
obj_db = cbind(obj_db,"out-dir"=create_out_dir(obj_db))                                                                               # add output directory

# check out-dir status (missing inputs, exists, up-to-date)
if (opt$verbose) write("Checking input/output object status...",stderr())
obj_db = cbind(obj_db,"status"=apply(obj_db,1,check_status))

# format input objects
if (opt$verbose) write("Formatting input objects...",stderr())
f = function(s) paste(paste("'",gsub(',',' ',unlist(strsplit(s,split=' '))),"'",sep=''),collapse=' ')
obj_db[,"inp-object"] = sapply(as.vector(obj_db[,"inp-object"]),f)                                                                    # input objects: convert commas to spaces and add single quotes

# print output object matrix
#if (opt$verbose) { write("OUTPUT OBJECT TABLE: ",stderr()); write.table(obj_db,quote=FALSE,col.names=TRUE,row.names=TRUE,sep='\t',stderr()) }

# generate shell commands
if (opt$verbose) write("Generating shell commands...",stderr())
obj_db = cbind(obj_db,"command"=generate_command(obj_db))                                                                             # add command line

# convert to data frame
obj_db = as.data.frame(obj_db)

# status and stats
compute = (obj_db[,"status"]=="new")
outofdate = (obj_db[,"status"]=="out-of-date")
if (opt$verbose) {
  write("Generating summary...",stderr())
  write(paste("- input objects = ",nrow(inp_db),sep=''),stderr())
  write(paste("- selected input objects = ",length(selected),sep=''),stderr())
  write(paste("- output objects = ",nrow(obj_db),sep=''),stderr())
  write(paste("- up-to-date output objects = ",sum(obj_db[,"status"]=="up-to-date"),sep=''),stderr())
  write(paste("- missing inputs = ",sum(obj_db[,"status"]=="missing-inputs"),sep=''),stderr())
  write(paste("- out-of-date output objects = ",sum(obj_db[,"status"]=="out-of-date"),sep=''),stderr())
  write(paste("- new output objects = ",sum(obj_db[,"status"]=="new"),sep=''),stderr())
  write(paste("- output objects to be computed = ",sum(compute),sep=''),stderr())
}

# arrange columns of object db in alphabetical order
obj_db = obj_db[,order(colnames(obj_db))]

# save data
if (opt$verbose) write("Saving database...",stderr())
save.image(file=paste(out_dir,"db.RData",sep="/"))
write.table(obj_db[compute,"command"],col.names=FALSE,row.names=FALSE,sep='\t',quote=FALSE,file=paste(out_dir,"run",sep="/"))
write.table(obj_db[outofdate,"command"],col.names=FALSE,row.names=FALSE,sep='\t',quote=FALSE,file=paste(out_dir,"run.outofdate",sep="/"))
write.table(obj_db,col.names=TRUE,row.names=FALSE,sep='\t',quote=FALSE,file=paste(out_dir,"db.tsv",sep="/"))

# done
if (opt$'verbose'==TRUE) write('Done.',stderr())
quit(save='no')




