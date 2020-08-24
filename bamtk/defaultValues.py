#########################################################################################
#                                                                                       #
# defaultValues.py - store default values used in many places in CheckM                 #
#                                                                                       #
#########################################################################################
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

class DefaultValues():
    """Default values for filenames and common constants."""

    FEATURES_ABUNDANCE_FILES = ['features_reads_raw_abundance.tsv','features_reads_normalised_abundance.tsv','features_reads_relative_abundance.tsv','features_base_raw_abundance.tsv','features_base_normalised_abundance.tsv','features_base_relative_abundance.tsv' ]

    ANNOTATE_ABUNDANCE_FILES = ['annotate_reads_raw_abundance.tsv', 'annotate_reads_normalised_abundance.tsv', 'annotate_reads_relative_abundance.tsv', 'annotate_base_raw_abundance.tsv', 'annotate_base_normalised_abundance.tsv' , 'annotate_base_relative_abundance.tsv']

