#########################################################################################
#                                                                                       #
#    This program is free software: you can redistribute it and/or modify               #
#    it under the terms of the GNU General Public License as published by               #
#    the Free Software Foundation, either version 3 of the License, or                  #
#    (at your option) any later version.                                                #
#                                                                                       #
#    This program is distributed in the hope that it will be useful,                    #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of                     #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                      #
#    GNU General Public License for more details.                                       #
#                                                                                       #
#    You should have received a copy of the GNU General Public License                  #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.               #
#                                                                                       #
#########################################################################################

import os
import sys

import argparse
import inspect

from biolib.common import make_sure_path_exists
from biolib.logger import logger_setup
from biolib.misc.custom_help_formatter import CustomHelpFormatter

from bamtk.main import OptionsParser

def print_help():
    """Help function"""

    print ('')
    print ('            ...::: BAM-Tk v' + version() + ' :::...')
    print ('''\

    Making Matrix: 
        mm_features             -> making bam features matrix 
        mm_annotated_features   -> making annoted features abundance matrix 
        mm_wf                   -> making bam features and annoted features abundance matrix 


Use: bamtk <command> -h for command specific help. 

Feature requests or bug reports can be sent to Corentin Hochart (corentin.hochart.pro@gmail.com)
or posted on GitHub (https://github.com/meb-team/Tools/).
    ''')


def version():
    import bamtk
    versionFile = open(os.path.join(bamtk.__path__[0], 'VERSION'))
    return versionFile.readline().strip()


def main():

    # initialize the option parser
    parser = argparse.ArgumentParser(add_help=False,
        description="BAM-Tk is a software toolkit for dealing with Binary Alignment Map (BAM) files.",
        epilog="Written by Corentin Hochart (corentin.hochart.pro@gmail.com), " +
        "UMR CNRSS 6023 Laboratoire Genome et Environement (LMGE), " +
        "as part of the [ANR Eureka](https://anr.fr/Projet-ANR-14-CE02-0004) project." +
        "Released under the terms of the GNU General Public License v3. " +
        "bamtk version %s." % version())
    subparsers = parser.add_subparsers(help="--", dest='subparser_name')

    # pathway reconstruction 
    mm_featuresparser = subparsers.add_parser('mm_features',
                                            description='') 
    mm_featuresparser.add_argument('faidx',help='samtools fasta index of the reference')
    mm_featuresparser.add_argument('bam_list',help='list of bam format alignement file(s) path')
    mm_featuresparser.add_argument('output_dir',help='directory to write output files')
    mm_featuresinput_argument = mm_featuresparser.add_argument_group('optional input arguments')
    mm_featuresinput_argument.add_argument('-x','--extension', help='bam file prefix',default='bam')
    mm_featuresinput_argument.add_argument('-fx','--faidx_extension', help='faidx file prefix',default='fasta.fai')
    mm_featuresinput_argument.add_argument('-t','--threads', help='threads number for "samtools view"',default='2')
    mm_featuresinput_argument.add_argument('-Q','--mapQ',help='only include reads with mapping quality >= INT [10]',default='10')
    mm_featuresinput_argument.add_argument('-i','--id_cutoff',help='only include reads with identity >= INT [0]',default=0)
    mm_featuresinput_argument.add_argument('-m','--merge',help='merge features abundance by field',action='store_true')
    mm_featuresinput_argument.add_argument('-s','--separator',help='filed separator for -m/--merge argument',default='.')
    mm_featuresinput_argument.add_argument('-g','--genome',help='sum abundance of all features',action='store_true')
    mm_featuresoutput_argument = mm_featuresparser.add_argument_group('optional output arguments')
    mm_featuresoutput_argument.add_argument('-n','--feature_normalisation',help="get the number of features per X reads [Default: 1000000]",default=1000000,type=int)
    mm_featuresoutput_argument.add_argument('-sn','--feature_size_normalisation',help="get the number of features per X bases [Default: 1000]",default=1000,type=int)
    mm_featuresoutput_argument.add_argument('-f','--discard_feature_length_normalisation',help="discard feature length normalisation for base count abundance output",action='store_true')
    mm_featuresoutput_argument.add_argument('-l','--discard_library_size_normalisation',help="discard library size normalisation for reads and bases count abundance output",action='store_true')
    mm_featuresoutput_argument.add_argument('-lsn','--library_size_normalisation',help="library size normalisation by total number of reads count or by number of aligned reads ",choices=['total','aligned'],default='total')
    mm_featuresoutput_argument.add_argument('--removed',help="removed features who do not appears in samples (sum of abundance through sample = 0)",action='store_true')
    mm_featuresparser.add_argument('--silent', help='suppress output of logger', action='store_true')
    mm_featuresparser.add_argument('--force_overwrite', help='force overwriting of output directory', action="store_true", default=False)
    mm_featuresparser.add_argument('--version',help='print version and exit',action='version',version='bamtk '+ version())


    mm_annotated_features_parser = subparsers.add_parser('mm_annotated_features',
                                            description='')
    mm_annotated_features_parser.add_argument('features_dir',help='directory specified during features command')
    mm_annotated_features_parser.add_argument('features_annotation',help='features annotation file in tabular format')
    mm_annotated_features_parser.add_argument('annotation_description',help='annotation description file in tabular format')
    mm_annotated_features_input_argument = mm_annotated_features_parser.add_argument_group('optional input arguments')
    mm_annotated_features_input_argument.add_argument('--library_size', help="Tabular file with sample library size to produce normalised count matrix")
    mm_annotated_features_output_argument = mm_annotated_features_parser.add_argument_group('optional output arguments')
    mm_annotated_features_output_argument.add_argument('-f','--feature_normalisation',help="get the number of features per X reads [Default: 1000000]",default=1000000,type=int)
    mm_annotated_features_output_argument.add_argument('--removed',help="removed features who do not appears in samples (sum of abundance through sample = 0)",action='store_true')                                            
    
    mm_annotated_features_parser.add_argument('--silent', help='suppress output of logger', action='store_true')
    mm_annotated_features_parser.add_argument('--force_overwrite', help='force overwriting of output directory', action="store_true", default=False)

    mm_wf_parser = subparsers.add_parser('mm_wf',
                                            description='Run features and annotate_features command',
                                            epilog='bamtk mm_wf ./file.fai ./bam_list.tsv ./features2annotation.tsv ./annotationDescription.tsv ./output')
    mm_wf_parser.add_argument('faidx',help='samtools fasta index of the reference')
    mm_wf_parser.add_argument('bam_list',help='list of bam format alignement file(s) path ')
    mm_wf_parser.add_argument('features_annotation',help='features annotation file in tabular format')
    mm_wf_parser.add_argument('annotation_description',help='annotation description file in tabular format')
    mm_wf_parser.add_argument('output_dir',help='directory to write output files')
    mm_wf_input_argument = mm_wf_parser.add_argument_group('optional input arguments')
    mm_wf_input_argument.add_argument('-x','--extension', help='bam file prefix',default='bam')
    mm_wf_input_argument.add_argument('-fx','--faidx_extension', help='faidx file prefix',default='fasta.fai')
    mm_wf_input_argument.add_argument('-t','--threads', help='threads number for "samtools view"',default='2')
    mm_wf_input_argument.add_argument('-Q','--mapQ',help='only include reads with mapping quality >= INT [10]',default='10')
    mm_wf_input_argument.add_argument('-i','--id_cutoff',help='only include reads with identity >= INT [0]',default=0)
    mm_wf_output_argument = mm_wf_parser.add_argument_group('optional output arguments')
    mm_wf_output_argument.add_argument('-f','--feature_normalisation',help="get the number of features per X reads [Default: 1000000]",default=1000000,type=int)
    mm_wf_output_argument.add_argument('-g','--discard_gene_length_normalisation',help="discard gene length normalisation for base count abundance output",action='store_true')
    mm_wf_output_argument.add_argument('--removed',help="removed features who do not appears in samples (sum of abundance through sample = 0)",action='store_true')
    mm_wf_parser.add_argument('--silent', help='suppress output of logger', action='store_true')
    mm_wf_parser.add_argument('--force_overwrite', help='force overwriting of output directory', action="store_true", default=False)
    mm_wf_parser.add_argument('--version',help='print version and exit',action='version',version='bamtk '+ version())

    # get and check options
    args = None
    if(len(sys.argv) == 1 or sys.argv[1] == '-h' or sys.argv == '--help'):
        print_help()
        sys.exit(0)
    else:
        args = parser.parse_args()

    try:
        logger_setup(args.output_dir, "bamtk.log", "bamtk", version(), args.silent)
    except:
        logger_setup(None, "bamtk.log", "bamtk", version(), args.silent)

    try:
        parser = OptionsParser()
        if(False):
            import cProfile
            cProfile.run('parser.parse_options(args)')
        else:
            parser.parse_options(args)
    except SystemExit:
        print ('Unrecoverable error.')
    except:
        print ("\nUnexpected error:", sys.exc_info()[0])
        raise

if __name__ == '__main__':
    main()