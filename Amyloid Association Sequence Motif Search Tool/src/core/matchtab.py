from Bio import ExPASy
from Bio import SeqIO
import re
import tkinter as tk
import tkinter.ttk as ttk
import requests as r
from io import StringIO


# Set up the tab showing a match
class MatchTab(tk.Frame):
    def __init__(self, master, sequence, match):
        tk.Frame.__init__(self, master)
        self.master = master
        self.sequence = sequence
        self.match = match
        self.setConfigs()
        self.createWidgets()

    def setConfigs(self):
        self.config(pady=2, padx=2)

    def createWidgets(self):
        #
        self.sequenceFrame = tk.Frame(self)
        self.sequenceFrame.grid(row=1, column=1, sticky="nesw")
        #
        sequenceTextScrollY = tk.Scrollbar(self.sequenceFrame)
        sequenceTextScrollY.pack(side=tk.RIGHT, fill=tk.Y)
        self.sequenceText = tk.Text(self.sequenceFrame, state="disabled", height=12,
                                    yscrollcommand=sequenceTextScrollY.set)
        self.sequenceText.pack(fill=tk.X)
        sequenceTextScrollY.config(command=self.sequenceText.yview)
        self.sequenceText.config(state="normal")
        self.sequenceText.insert(tk.INSERT, self.sequence.seq)
        #
        # Highlight the mathc in the sequence
        if "highlight" in self.sequenceText.tag_names():
            self.sequenceText.tag_delete("highlight")
        i = len(self.match.group(1))
        index = "1.0"
        while True:
            index = self.sequenceText.search(self.match.group(1), index, nocase=1, stopindex="end")
            if index:
                index2 = self.sequenceText.index("%s+%dc" % (index, i))
                self.sequenceText.tag_add("highlight", index, index2)
                self.sequenceText.tag_config("highlight", background="yellow")
                index = index2
            else:
                break
        self.sequenceText.config(state="disabled")
        #
        self.infoFrame = tk.Frame(self)
        self.infoFrame.grid(row=2,column=1,pady=(5,0),sticky="nsew")
        #
        self.matchFrame = tk.Frame(self.infoFrame)
        self.matchFrame.pack()
        #
        self.matchLabel = tk.Label(self.matchFrame, text="Match:")
        self.matchLabel.grid(row=1,column=1)
        #
        self.matchText = tk.Text(self.matchFrame,state="normal", height=1, width=20)
        self.matchText.insert(tk.INSERT, self.match.group(1))
        self.matchText.config(state="disabled")
        self.matchText.grid(row=2,column=1)
        #
        self.posFrame = tk.Frame(self.infoFrame)
        self.posFrame.pack()
        #
        self.posLabel = tk.Label(self.posFrame, text="At position:")
        self.posLabel.grid(row=1,column=1,columnspan=3)
        #
        self.posStartEntry = tk.Text(self.posFrame,state="normal", height=1, width=20)
        self.posStartEntry.insert(tk.INSERT, self.match.start(1)+1)
        self.posStartEntry.config(state="disabled")
        self.posStartEntry.grid(row=2,column=1)
        #
        self.posToLabel = tk.Label(self.posFrame, text="to")
        self.posToLabel.grid(row=2,column=2)
        #
        self.posEndEntry = tk.Text(self.posFrame,state="normal", height=1,width=20)
        self.posEndEntry.insert(tk.INSERT, self.match.end(1))
        self.posEndEntry.config(state="disabled")
        self.posEndEntry.grid(row=2,column=3)
        
