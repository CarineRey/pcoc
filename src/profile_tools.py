#  profile_tools.py
#
#  Copyright 2019 Carine Rey <carine.rey@ens-lyon.fr>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#

import os, sys
import argparse

import logging
logger = logging.getLogger("pcoc.profile_tools")

import pandas as pd


class Profiles(object):
    def __init__(self,input_profiles):
        if input_profiles in ["C10","C60"]:
            self.name = input_profiles
            self.filename = None
            self.nb_cat = int(input_profiles.replace("C",""))
            self.formatted_frequencies_filename = None
            self.df = None
        elif os.path.isfile(input_profiles):
            self.name = "from_file (%s)" %(input_profiles)
            self.filename = input_profiles
            self.nb_cat = 0
            self.formatted_frequencies_filename = None
            self.df = None
        else:
            logger.error(" invalid choice: %s (choose from 10, 60 or a csv filename containing amino acid frequencies for each profil)", value)
            sys.exit(1)

    def format_profile_file(self):
        try:
            self.df = pd.read_csv(self.filename, index_col = 0)
            # Check row names:
            index = list(self.df.index)
            aa_ordered = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
            if len(index) == 20 and set(index) == set(aa_ordered):
                self.df.reindex(aa_ordered)
            else:
                self.df=pd.DataFrame()
            self.nb_cat = self.df.shape[1]
        except Exception as exc:
            logger.error(str(exc))
            self.df=pd.DataFrame()
        if self.df.empty:
            logger.error("%s is not well formated, it must be a csv file and row names must be A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V.", self.filename)
            sys.exit(1)



def check_profiles(x_profiles):
    profiles = Profiles(x_profiles)
    if profiles.name == "from_file":
        profiles.format_profile_file()
    return profiles

