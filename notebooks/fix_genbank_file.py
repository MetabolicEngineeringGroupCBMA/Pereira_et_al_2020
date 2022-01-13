#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pydna.genbankfixer import gbtext_clean

with open("pTDH3-tc3-6xHA.gb") as f:
    gb = f.read()

cleaned = gbtext_clean(gb)

from pydna.readers import read

seq = read(cleaned.gbtext)

seq.write("pTDH3-tc3-6xHA_fixed.gb")
