#####################
# BASIC
###################

umask 022

export LC_ALL=C

#setenv MANPATH /usr/local/man:/usr/man
#setenv _POSIX2_VERSION 199209
#setenv LC_COLLATE "C"

#set noclobber
#set autolist
#set autoexpand
#set nohup
#set nobeep
#unset autocorrect
#unset correct


#####################
# PATH
###################

export PATH=$PATH:/ifs/home/at570/disk1/Resources/Bin/linux


#####################
# PROMPT
###################

export PS1="[\u@\h:\w ] $ "


#####################
# ALIASES
###################

# Aliases for date/time
alias date='date "+%A %d %B %Y, %T"'

# Aliases for ls
alias ll='ls -l --color=auto' 
alias lla='ls -aglF --color=auto'

# Aliases for mv, cp, rm
alias mv='mv -i'
alias cp='cp -i'
alias rm='rm -i'

# Miscelleaneous
alias c=clear
alias j=jobs
alias xterm="xterm -fg black -bg white -sb"
alias joint="join -t'	'"

# ssh/sftp
alias ssh-slinky="ssh -X tsirigos@slinky.cs.nyu.edu"
alias ssh-phoenix="ssh -X at570@phoenix.med.nyu.edu"
alias ssh-aifantis="ssh -X at570@10.193.37.34"
alias ssh-nyuhpc="ssh -X at570@hpc.nyu.edu"
alias ftp-slinky="sftp tsirigos@slinky.cs.nyu.edu"


#####################
# MISC
###################

# Added by the mirnylab install script.
export PYTHONPATH="$PYTHONPATH:$HOME/disk1/Resources/downloads/mirnylib"

# Added by the mirnylab install script.
export PYTHONPATH="$PYTHONPATH:$HOME/disk1/Resources/downloads/hiclib/src"



