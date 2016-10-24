# Set the path to the pipelines folder.
# Modify as needed.
# $TOOLS is an environment variable set to the tools folder in my .bash_profile.
# $TOOLs can be replacaed by the full path to the tools folder.
export PIPELINES_FOLDER=$TOOLS/abcngspipelines

#######################################################
# Sets paths to Python scripts used in NGS pipelines. #
#######################################################

PATH=$PIPELINES_FOLDER/alignment:$PATH
PATH=$PIPELINES_FOLDER/bisseq:$PATH
PATH=$PIPELINES_FOLDER/bischipseq:$PATH
PATH=$PIPELINES_FOLDER/chipseq:$PATH
PATH=$PIPELINES_FOLDER/chipseq/homer:$PATH
PATH=$PIPELINES_FOLDER/exomeseq:$PATH
PATH=$PIPELINES_FOLDER/nanuq:$PATH
PATH=$PIPELINES_FOLDER/pipeline:$PATH
PATH=$PIPELINES_FOLDER/quality_control:$PATH
PATH=$PIPELINES_FOLDER/report:$PATH
PATH=$PIPELINES_FOLDER/rnaseq:$PATH
PATH=$PIPELINES_FOLDER/rnaseq/deseq:$PATH
PATH=$PIPELINES_FOLDER/samtools:$PATH
PATH=$PIPELINES_FOLDER/smallrnaseq:$PATH
PATH=$PIPELINES_FOLDER/sra:$PATH
PATH=$PIPELINES_FOLDER/trimming:$PATH
PATH=$PIPELINES_FOLDER/utils:$PATH
PATH=.:$PATH
export PATH

PYTHONPATH=$PIPELINES_FOLDER/utils:$PYTHONPATH
export PYTHONPATH
