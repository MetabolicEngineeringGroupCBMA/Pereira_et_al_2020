# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.1
#   kernelspec:
#     display_name: Python [conda env:bjorn39]
#     language: python
#     name: conda-env-bjorn39-py
# ---

# ## Jupyter notebook
#
# comment

from pydna.amplify import Anneal
from pydna.amplify import pcr
from pydna.assembly import Assembly
from pydna.genbank import genbank
from pydna.genbank import Genbank
from pydna.download import download_text
from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from pydna.readers import read
from pydna.readers import read_primer
from pydna.parsers import parse
from pydna.parsers import parse_primers
from pydna.editor import ape
from pydna.design import primer_design
from pydna.design import assembly_fragments
from pydna.design import circular_assembly_fragments
from pydna.utils import eq
from pydna.genbankfixer import gbtext_clean

pl = parse_primers("""

>671_sc_acc1-Tc2: (83-mer)
CGATACGATACGACACGATACGATACGACACGCTACTATAGCATAGGCCACTAGTGGATCTG

>670_sc_acc1-Tc1B: (82-mer)
CCATCTTCTGTGGAGAAGACTCGAATAAGCTTTCTTCGCTCATATGTTCTCGAGGCCTAGG

>577_crp585-557 (29-mer)
gttctgatcctcgagcatcttaagaattc

>567_pCAPsAjiIF (23-mer)
GTCggctgcaggtcactagtgag

>468_pCAPs_release_fw (25-mer) 79.66 same as 560
gtcgaggaacgccaggttgcccact

>1259_ScACC1middleRV
CCTTCGTGAACTCTAATATCTCC

>1260_ScACC1middleFW
GCTCAAGTCTATATTCGTCG

>467_pCAPs_release_re (31-mer)
ATTTAAatcctgatgcgtttgtctgcacaga

>568_pCAPsAjiIR (22-mer)
GTGCcatctgtgcagacaaacg

>578_crp42-70 (29-mer)
gttcttgtctcattgccacattcataagt

>577_crp585-557 (29-mer)
gttctgatcctcgagcatcttaagaattc

>567_pCAPsAjiIF (23-mer)
GTCggctgcaggtcactagtgag

>468_pCAPs_release_fw (25-mer) 79.66 same as 560
gtcgaggaacgccaggttgcccact

>779_YlACC1_3445_rv (25-mer)
ACAAAGCAGACGACATGGTAGGCAG

>780_YlACC1_3305_fwd (25-mer)
TCTTTGCCCACGATGATCCCTGGAT

>467_pCAPs_release_re (31-mer)
ATTTAAatcctgatgcgtttgtctgcacaga

>568_pCAPsAjiIR (22-mer)
GTGCcatctgtgcagacaaacg

>578_crp42-70 (29-mer)
gttcttgtctcattgccacattcataagt

>577_crp585-557 (29-mer)
gttctgatcctcgagcatcttaagaattc

>567_pCAPsAjiIF (23-mer)
GTCggctgcaggtcactagtgag

>1258_ScACC1_fw
cccactttctcactagtgacctgcagccgacaaATGAGCGAAGAAAGCT

>1257_ScACC1_rv
aaatcctgatgcgtttgtctgcacagatggcacTTATTTCAAAGTCTTCAACAAT

>568_pCAPsAjiIR (22-mer)
GTGCcatctgtgcagacaaacg

>578_crp42-70 (29-mer)
gttcttgtctcattgccacattcataagt

>468_pCAPs_release_fw (25-mer) 79.66 same as 560
gtcgaggaacgccaggttgcccact

>1259_ScACC1middleRV
CCTTCGTGAACTCTAATATCTCC

>1260_ScACC1middleFW
GCTCAAGTCTATATTCGTCG

>467_pCAPs_release_re (31-mer)
ATTTAAatcctgatgcgtttgtctgcacaga

>577_crp585-557 (29-mer)
gttctgatcctcgagcatcttaagaattc

>1259_ScACC1middleRV
CCTTCGTGAACTCTAATATCTCC

>1260_ScACC1middleFW
GCTCAAGTCTATATTCGTCG

>623_ScTDH3tpr_PacI (33-mer)
taattaaTTTGTTTGTTTATGTGTGTTTATTCG

>1123_New775
gtgcaatgcggccgctGAC

>779_YlACC1_3445_rv (25-mer)
ACAAAGCAGACGACATGGTAGGCAG

>780_YlACC1_3305_fwd (25-mer)
TCTTTGCCCACGATGATCCCTGGAT

>578_crp42-70 (29-mer)
gttcttgtctcattgccacattcataagt

>577_crp585-557 (29-mer)
gttctgatcctcgagcatcttaagaattc

>567_pCAPsAjiIF (23-mer)
GTCggctgcaggtcactagtgag

>781_YlACC1f_YPK (63-mer)
GCCAGGTTGCCCACTTTCTCACTAGTGACCTGCAGCCCACatgcgactgcaattgaggac
act

>505_YlACC1r_CpoI
TCCGtcacaaccccttgagcagctca

>504_YlACC1f_SgsI
CCAAatgcgactgcaattgaggacact

>782_YlACC1r_YPK (62-mer)
TAAATCCGGATATCCTGATGCGTTTGTCTGCACAGATGACtcacaaccccttgagcagct
ca

>568_pCAPsAjiIR (22-mer)
GTGCcatctgtgcagacaaacg

>578_crp42-70 (29-mer)
gttcttgtctcattgccacattcataagt

>577_crp585-557 (29-mer)
gttctgatcctcgagcatcttaagaattc

>779_YlACC1_3445_rv (25-mer)
ACAAAGCAGACGACATGGTAGGCAG

>780_YlACC1_3305_fwd (25-mer)
TCTTTGCCCACGATGATCCCTGGAT

>623_ScTDH3tpr_PacI (33-mer)
taattaaTTTGTTTGTTTATGTGTGTTTATTCG

>1123_New775
gtgcaatgcggccgctGAC

>1259_ScACC1middleRV
CCTTCGTGAACTCTAATATCTCC

>1260_ScACC1middleFW
GCTCAAGTCTATATTCGTCG

>578_crp42-70 (29-mer)
gttcttgtctcattgccacattcataagt

""")


















"""
415
467
468
504
505
564
567
568
577
578
586
622
623
670
671
698
779
780
781
782
1123
1257
1258
1259
1260
1282
"""


others = [415,564,586,622,698,1282]

from pydna.myprimers import PrimerList

plst = PrimerList()

plst.pydna_code_from_indices(others)



"""

>415_ScTDH3tpf (29-mer)
TTAAATAATAAAAAACACGCTTTTTCAGT

>564_YlACC1_628_R (25-mer)
tcgtccactccggttccagaccacg

>586_YlACC1_6264_F (20-mer)
CTCCTCTCTCAAGAAGCAGC

>622_ScPGI1tpr_PacI (28-mer)
taattaaTTTTAGGCTGGTATCTTGATT

>698_sc_acc1-B1: (21-mer)
ACCTGGCACTTCAATGTATTG

>1282_sc_acc1-T
GCGACCATGACAATGCTATTGATGG

"""

467
468
504
505
564
567
568
577
578
586
622
623
670
671
698
779
780
781
782
1123
1258
1259
1260
1282




"""
415
467	467_pCAPs_release_re	ATTTAAatcctgatgcgtttgtctgcacaga
468	468_pCAPs_release_fw	gtcgaggaacgccaggttgcccact
504	504_YlACC1f_SgsI	CCAAatgcgactgcaattgaggacact
505	505_YlACC1r_CpoI	TCCGtcacaaccccttgagcagctca
564
567	567_pCAPsAjiIF	GTCggctgcaggtcactagtgag
568	568_pCAPsAjiIR	GTGCcatctgtgcagacaaacg
577	577_crp585-557	gttctgatcctcgagcatcttaagaattc
578	578_crp42-70	gttcttgtctcattgccacattcataagt
586
622
623	623_ScTDH3tpr_PacI	taattaaTTTGTTTGTTTATGTGTGTTTATTCG
670	670_sc_acc1-Tc1B:	CCATCTTCTGTGGAGAAGACTCGAATAAGCTTTCTTCGCTCATATGTTCTCGAGGCCTAGG
671	671_sc_acc1-Tc2:	CGATACGATACGACACGATACGATACGACACGCTACTATAGCATAGGCCACTAGTGGATCTG
698
779	779_YlACC1_3445_rv	ACAAAGCAGACGACATGGTAGGCAG
780	780_YlACC1_3305_fwd	TCTTTGCCCACGATGATCCCTGGAT
781	781_YlACC1f_YPK	GCCAGGTTGCCCACTTTCTCACTAGTGACCTGCAGCCCACatgcgactgcaattgaggacact
782	782_YlACC1r_YPK	TAAATCCGGATATCCTGATGCGTTTGTCTGCACAGATGACtcacaaccccttgagcagctca
1123	1123_New775	gtgcaatgcggccgctGAC
1257	1257_ScACC1_rv	aaatcctgatgcgtttgtctgcacagatggcacTTATTTCAAAGTCTTCAACAAT
1258	1258_ScACC1_fw	cccactttctcactagtgacctgcagccgacaaATGAGCGAAGAAAGCT
1259	1259_ScACC1middleRV	CCTTCGTGAACTCTAATATCTCC
1260	1260_ScACC1middleFW	GCTCAAGTCTATATTCGTCG
1282
"""




