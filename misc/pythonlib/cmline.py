from Bio import Application 
from Bio.Application import _Option 
class NeedleCommandline(Application.AbstractCommandline):
    """Commandline object for the needle program from EMBOSS.
    """
    def __init__(self, cmd = "needle", **kwargs):

        self.program_name = cmd
        self.parameters = \
            [_Option(["-asequence", "-a"], ["input", "file"], None, 1,
                     "First sequence to align"),
             _Option(["-bsequence", "-b"], ["input", "file"], None, 1,
                     "Second sequence to align"),
             _Option(["-gapopen", "-go"], ["input"], None, 1,
                     "Gap open penalty"),
             _Option(["-gapextend", "-ge"], ["input"], None, 1,
                     "Gap extension penalty"),
             _Option(["-outfile", "-o"], ["output", "file"], None, 1,
                     "Output file for the alignment"),
             _Option(["-sreverse1", "-rev1"], ["input"], None, 0,
                     "Reverse complement asequence"),
             _Option(["-sreverse2", "-rev2"], ["input"], None, 0,
                     "Reverse complement asequence"),
             _Option(["-datafile"], ["input", "file"], None, 0,
                     "Matrix file"),
             _Option(["-similarity"], ["input"], None, 0,
                     "Display percent identity and similarity"),
             _Option(["-nosimilarity"], ["input"], None, 0,
                     "Do not display percent identity and similarity"),
             _Option(["-aformat"], ["input"], None, 0,
                     "Display output in a different specified output format"),
             _Option(["-auto"], ["input"], None, 0,
                     "Turns off prompt"),
             _Option(["-filter"], ["input"], None, 0,
                     "Reads asequence from stdin"),
             _Option(["-ausashow3", "-usa"], ["input"], None, 0,
                     "Show the full USA in the alignment"),
             _Option(["-adesshow3", "-des"], ["input"], None, 0,
                     "Show description in the header")
             ]
        
