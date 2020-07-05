from Bio import ExPASy
from Bio import SeqIO
import re
import tkinter as tk
import tkinter.ttk as ttk
import requests as r
from io import StringIO

from matchtab import MatchTab


# A windows that displays the matches for a protein. Each match is in a tab of its own.
# The full sequence of the protein is shown for each watch with the match itself highlighted in yellow.
# Start and end positions of the match is also shown alongside the match itself.
class InfoViewer(tk.Frame):
    def __init__(self, master, code, sequence, matches, num):
        tk.Frame.__init__(self, master)
        self.master = master
        self.code = code
        self.sequence = sequence[0]
        self.matches = matches
        self.num = num
        self.setConfigs()
        self.createWidgets()
        self.master.config(pady=5, padx=5)

    # Set configs of the window. Title displays the line number selected in the main window and the UNIPROT code of the protein used.
    def setConfigs(self):
        self.master.title("[%s] %s"% (self.num, self.code))
        self.master.resizable(False,False)

    # Set up widgets for window. Displays UNIPROT code and creats tabs for each matches.
    def createWidgets(self):
        self.titleLabel = tk.Label(self, height=2, text="Displaying motif matches for Uniprot Accession Code:\n" + self.code)
        self.titleLabel.grid(row=1, column=1, sticky="ew")
        self.matchNote = ttk.Notebook(self)
        self.matchNote.grid(row=2,column=1, sticky="nwes")
        num = 0
        for match in self.matches:
            num += 1
            self.matchNote.add(MatchTab(self.matchNote, self.sequence, match), text="Match " + str(num))
